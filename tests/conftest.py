import pytest
from testfixtures import TempDirectory


@pytest.fixture()
def create_bad_xlsx():
    # setup temp folder for fake .xlsx file
    temp = TempDirectory()
    temp_dir = temp.path

    # create dummy .xlsx
    # build dictionaries, use pandas to write .xlsx to temp
    # yield temp folder

    # breakdown temp folder
    #   delete temp and .xlsx

@pytest.fixture()
def create_good_xlsx():
    pass