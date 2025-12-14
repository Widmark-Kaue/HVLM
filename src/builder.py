from pathlib import Path
from configparser import ConfigParser
from avlwrapper import Configuration

################################################################
# AVL and Project Path
################################################################

# avlPath = Path('bin','avl352.exe').absolute()
avlPath = Path('bin','avl.exe').absolute()
if not avlPath.exists():
    # avlPath = Path('..', 'bin', 'avl352.exe').absolute()
    avlPath = Path('..', 'bin', 'avl.exe').absolute()

PROJECT = avlPath.parent.parent

cfg_path = PROJECT.joinpath('config.cfg')

################################################################
# Setup to run avl
################################################################
def config_file(PATH_AVL):
    avl_path = PATH_AVL

    config = ConfigParser()

    config['environment'] = {
        'Executable': str(avl_path.absolute()),
        'PrintOutput': 'no',
        'GhostscriptExecutable': 'gs',
        'loglevel': 'INFO'
    }

    config['output'] = {
        'Totals': 'yes',
        'SurfaceForces': 'yes',
        'StripForces': 'yes',
        'ElementForces': 'yes',
        'BodyAxisDerivatives': 'yes',
        'StabilityDerivatives': 'yes',
        'HingeMoments': 'yes',
        'StripShearMoments': 'yes',
    }

    with open(cfg_path, 'w') as config_file:
        config.write(config_file)

config_file(avlPath)
my_config = Configuration(cfg_path)
