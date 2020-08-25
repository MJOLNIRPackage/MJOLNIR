def pytest_addoption(parser):
    parser.addoption("--quick", action="store_true")


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    option_value = metafunc.config.option.quick
    if 'quick' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("quick", [option_value])
