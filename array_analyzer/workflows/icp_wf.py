# bchhun, {2020-04-02}

from datetime import datetime
import glob
import time
import skimage.io as io
import pandas as pd

import array_analyzer.extract.image_parser as image_parser
import array_analyzer.extract.img_processing as img_processing
import array_analyzer.extract.txt_parser as txt_parser
from array_analyzer.load.debug_images import *
from array_analyzer.transform.property_filters import *
import array_analyzer.transform.point_registration as registration

FIDUCIALS = [(0, 0), (0, 1), (0, 5), (7, 0), (7, 5)]
FIDUCIALS_IDX = [0, 5, 6, 30, 35]
SCENION_SPOT_DIST = 82


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
        # finding center of well and cropping
        cx, cy, r, well_mask = image_parser.find_well_border(
            image,
            segmethod='otsu',
            detmethod='region',
        )
        im_crop = image_parser.crop_image(image, cx, cy, r, border_=0)
        # Estimate and remove background
        background = img_processing.get_background(im_crop, fit_order=2)
        im_crop = (im_crop / background * np.mean(background)).astype(np.uint8)

        nbr_grid_rows, nbr_grid_cols = props_array.shape
        spot_coords = img_processing.get_spot_coords(
            im_crop,
            min_area=1000,
            min_thresh=25,
        )
        # Initial estimate of spot center
        mean_point, spot_dist = img_processing.grid_estimation(
            im=im_crop,
            spot_coords=spot_coords,
        )
        if spot_dist is None:
            spot_dist = SCENION_SPOT_DIST
        grid_coords = registration.create_reference_grid(
            mean_point=mean_point,
            nbr_grid_rows=nbr_grid_rows,
            nbr_grid_cols=nbr_grid_cols,
            spot_dist=SCENION_SPOT_DIST,
        )
        # Register grid to spots using iterative closest point
        t_matrix = registration.icp(
            source=grid_coords[FIDUCIALS_IDX],
            target=spot_coords,
        )
        if t_matrix is None:
            # Registration failed, try mean instead
            mean_point = tuple(np.mean(spot_coords, axis=0))
            grid_coords = registration.create_reference_grid(
                mean_point=mean_point,
                nbr_grid_rows=nbr_grid_rows,
                nbr_grid_cols=nbr_grid_cols,
                spot_dist=SCENION_SPOT_DIST,
            )
            t_matrix = registration.icp(
                source=grid_coords[FIDUCIALS_IDX],
                target=spot_coords,
            )
        if t_matrix is not None:
            # Transform coordinates
            reg_coords = np.squeeze(cv.transform(
                np.array([grid_coords]),
                t_matrix),
            )
        else:
            print("ICP failed, using initial estimate")
            reg_coords = grid_coords

        # od_well, i_well, bg_well = compute_od(props_array_placed, bgprops_array)
        #
        # pd_OD = pd.DataFrame(od_well)
        # pd_OD.to_excel(xlwriterOD, sheet_name=image_name[:-4])

        print("Time to register grid to {}: {:.3f} s".format(
            image_name,
            time.time() - start_time),
        )

        # SAVE FOR DEBUGGING
        if debug:
            # # Save image with spots
            im_roi = im_crop.copy()
            im_roi = cv.cvtColor(im_roi, cv.COLOR_GRAY2RGB)
            plt.imshow(im_roi)
            plt.plot(spot_coords[:,0],spot_coords[:,1],'rx',ms=12)
            plt.plot(grid_coords[:,0],grid_coords[:,1],'b+',ms=12)
            plt.plot(reg_coords[:,0],reg_coords[:,1],'g.',ms=10)

            write_name = image_name[:-4] + '_icp.jpg'
            figICP = plt.gcf()
            figICP.savefig(os.path.join(run_path, write_name))
            plt.show()
            # cv.imwrite flips the color identity. Confusing to write the diagnostic plot and interpret.
            # cv.imwrite(os.path.join(run_path, write_name), cv.cvtColor(im_roi, cv.COLOR_RGB2BGR))


            # well_path = os.path.join(run_path)
            # os.makedirs(well_path, exist_ok=True)
            # output_name = os.path.join(well_path, image_name[:-4])
            # im_bg_overlay = np.stack([background, im_crop, background], axis=2)
            #
            # # This plot shows which spots have been assigned what index.
            # plot_spot_assignment(od_well, i_well, bg_well,
            #                      im_crop, props_placed_by_loc, bgprops_by_loc,
            #                      image_name, output_name, params)
            #
            # #   save spots
            # save_all_wells(props_array, spot_ids, well_path, image_name[:-4])
            #
            # #   save a composite of all spots, where spots are from source or from region prop
            # save_composite_spots(im_crop, props_array_placed, well_path, image_name[:-4], from_source=True)
            # save_composite_spots(im_crop, props_array_placed, well_path, image_name[:-4], from_source=False)
            #
            # stop2 = time.time()
            # print(f"\ttime to save debug={stop2-stop}")
