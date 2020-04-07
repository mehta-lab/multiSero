# bchhun, {2020-04-02}

from datetime import datetime
import glob
import time
import skimage.io as io
import pandas as pd
import re
import cv2 as cv

import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.extract.txt_parser as txt_parser
from array_analyzer.load.debug_images import *
import array_analyzer.transform.point_registration as registration
from array_analyzer.transform.array_generation import build_centroid_binary_blocks

FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
FIDUCIALS_IDX = [0, 5, 6, 30, 35]
FIDUCIALS_IDX_8COLS = [0, 7, 8, 40, 47]
SCENION_SPOT_DIST = 82
# The expected standard deviations could be estimated from training data
# Just winging it for now
STDS = np.array([100, 100, .1, .001])  # x, y, angle, scale


def point_registration(input_folder_, output_folder_, debug=False):

    xml_path = glob.glob(input_folder_ + '/*.xml')
    if len(xml_path) > 1 or not xml_path:
        raise IOError("Did not find unique xml")
    xml_path = xml_path[0]

    # parsing .xml
    fiduc, spots, repl, params = txt_parser.create_xml_dict(xml_path)

    # creating our arrays
    spot_ids = txt_parser.create_array(params['rows'], params['columns'])
    antigen_array = txt_parser.create_array(params['rows'], params['columns'])

    # adding .xml info to these arrays
    spot_ids = txt_parser.populate_array_id(spot_ids, spots)

    antigen_array = txt_parser.populate_array_antigen(antigen_array, spot_ids, repl)

    # save a sub path for this processing run
    run_path = os.path.join(
        output_folder_,
        '_'.join([str(datetime.now().month),
                  str(datetime.now().day),
                  str(datetime.now().hour),
                  str(datetime.now().minute),
                  str(datetime.now().second)]),
    )
    xlwriterOD = pd.ExcelWriter(os.path.join(run_path, 'ODs.xlsx'))
    pdantigen = pd.DataFrame(antigen_array)
    pdantigen.to_excel(xlwriterOD, sheet_name='antigens')
    if debug:
        xlwriter_int = pd.ExcelWriter(os.path.join(run_path, 'intensities.xlsx'))
        xlwriter_bg = pd.ExcelWriter(os.path.join(run_path, 'backgrounds.xlsx'))

    if not os.path.isdir(run_path):
        os.mkdir(run_path)

    # ================
    # loop over images => good place for multiproc?  careful with columns in report
    # ================
    images = [file for file in os.listdir(input_folder_) if '.png' in file or '.tif' in file or '.jpg' in file]

    # remove any images that are not images of wells.
    wellimages = [file for file in images if re.match(r'[A-P][0-9]{1,2}', file)]

    # sort by letter, then by number (with '10' coming AFTER '9')
    wellimages.sort(key=lambda x: (x[0], int(x[1:-4])))

    for image_name in wellimages:
        start_time = time.time()
        image = image_parser.read_gray_im(os.path.join(input_folder_, image_name))

        props_array = txt_parser.create_array(
            params['rows'],
            params['columns'],
            dtype=object,
        )
        bgprops_array = txt_parser.create_array(
            params['rows'],
            params['columns'],
            dtype=object,
        )

        nbr_grid_rows, nbr_grid_cols = props_array.shape
        fiducials_idx = FIDUCIALS_IDX
        if nbr_grid_cols == 8:
            fiducials_idx = FIDUCIALS_IDX_8COLS

        spot_coords = image_parser.get_spot_coords(
            image,
            min_area=250,
            min_thresh=0,
        )

        # Initial estimate of spot center
        mean_point = tuple(np.mean(spot_coords, axis=0))
        grid_coords = image_parser.create_reference_grid(
            mean_point=mean_point,
            nbr_grid_rows=nbr_grid_rows,
            nbr_grid_cols=nbr_grid_cols,
            spot_dist=SCENION_SPOT_DIST,
        )
        fiducial_coords = grid_coords[fiducials_idx, :]

        particles = registration.create_gaussian_particles(
            mean_point=(0, 0),
            stds=STDS,
            scale_mean=1.,
            angle_mean=0.,
            nbr_particles=1000,
        )
        # Optimize estimated coordinates with iterative closest point
        t_matrix = registration.particle_filter(
            fiducial_coords=fiducial_coords,
            spot_coords=spot_coords,
            particles=particles,
            stds=STDS,
        )
        # Transform grid coordinates
        reg_coords = np.squeeze(cv.transform(np.array([grid_coords]), t_matrix))

        # Crop image
        im_crop, crop_coords = img_processing.crop_image_from_coords(
            im=image,
            grid_coords=reg_coords,
        )

        # Estimate and remove background
        background = img_processing.get_background(im_crop, fit_order=2)
        im_crop = (im_crop / background * np.mean(background)).astype(np.uint8)

        placed_spotmask = build_centroid_binary_blocks(
            reg_coords,
            im_crop,
            params,
        )
        spot_props = image_parser.generate_props(
            placed_spotmask,
            intensity_image_=im_crop,
        )
        bg_props = image_parser.generate_props(
            placed_spotmask,
            intensity_image_=background,
        )

        # unnecessary?  both receive the same spotmask
        spot_labels = [p.label for p in spot_props]
        bg_props = image_parser.select_props(
            bg_props,
            attribute="label",
            condition="is_in",
            condition_value=spot_labels,
        )
        props_placed_by_loc = image_parser.generate_props_dict(
            spot_props,
            params['rows'],
            params['columns'],
            min_area=100,
        )
        bgprops_by_loc = image_parser.generate_props_dict(
            bg_props,
            params['rows'],
            params['columns'],
            min_area=100,
        )

        props_array_placed = image_parser.assign_props_to_array(
            props_array,
            props_placed_by_loc,
        )
        bgprops_array = image_parser.assign_props_to_array(
            bgprops_array,
            bgprops_by_loc,
        )
        od_well, int_well, bg_well = image_parser.compute_od(
            props_array_placed,
            bgprops_array,
        )

        pd_OD = pd.DataFrame(od_well)
        pd_OD.to_excel(xlwriterOD, sheet_name=image_name[:-4])

        print("Time to register grid to {}: {:.3f} s".format(
            image_name,
            time.time() - start_time),
        )

        # ==================================

        # SAVE FOR DEBUGGING
        if debug:
            srt = time.time()
            well_path = os.path.join(run_path)
            os.makedirs(run_path, exist_ok=True)
            output_name = os.path.join(well_path, image_name[:-4])

            # Save spot and background intensities.
            pd_int = pd.DataFrame(int_well)
            pd_int.to_excel(xlwriter_int, sheet_name=image_name[:-4])
            pd_bg = pd.DataFrame(bg_well)
            pd_bg.to_excel(xlwriter_bg, sheet_name=image_name[:-4])

            # Save mask of the well, cropped grayscale image, cropped spot segmentation.
            # io.imsave(output_name + "_well_mask.png",
            #           (255 * well_mask).astype('uint8'))
            io.imsave(output_name + "_crop.png",
                      (255 * im_crop).astype('uint8'))
            # io.imsave(output_name + "_crop_binary.png",
            #           (255 * well_mask).astype('uint8'))

            # Evaluate accuracy of background estimation with green (image), magenta (background) overlay.
            im_bg_overlay = np.stack([background, im_crop, background], axis=2)
            io.imsave(output_name + "_crop_bg_overlay.png",
                      (255 * im_bg_overlay).astype('uint8'))

            # This plot shows which spots have been assigned what index.
            plot_spot_assignment(
                od_well, int_well,
                bg_well,
                im_crop,
                props_placed_by_loc,
                bgprops_by_loc,
                image_name,
                output_name,
                params,
            )
            # Save a composite of all spots, where spots are from source or from region prop
            save_composite_spots(
                im_crop,
                props_array_placed,
                well_path,
                image_name[:-4],
                from_source=True,
            )
            print(f"Time to save debug images: {time.time()-srt} s")

            # # Save image with spots
            im_roi = im_crop.copy()
            im_roi = cv.cvtColor(im_roi, cv.COLOR_GRAY2RGB)
            plt.imshow(im_roi)
            plt.plot(spot_coords[:, 0], spot_coords[:, 1], 'rx', ms=12)
            plt.plot(grid_coords[:, 0], grid_coords[:, 1], 'b+', ms=12)
            plt.plot(reg_coords[:, 0], reg_coords[:, 1], 'g.', ms=10)
            write_name = image_name[:-4] + '_registration.jpg'
            figICP = plt.gcf()
            figICP.savefig(os.path.join(run_path, write_name))
            plt.close(figICP)
            # cv.imwrite flips the color identity. Confusing to write the diagnostic plot and interpret.
            # cv.imwrite(os.path.join(run_path, write_name), cv.cvtColor(im_roi, cv.COLOR_RGB2BGR))
    if debug:
        xlwriter_int.close()
        xlwriter_bg.close()
    xlwriterOD.close()
