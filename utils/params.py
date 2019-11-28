import os
import re
from importlib.util import spec_from_file_location, module_from_spec


class Params:
    def configure(self, path):
        m = import_from_file(path)

        reg = re.compile(r"^__.+__$")  # Matches magic methods
        for name, value in m.__dict__.items():
            if reg.match(name):
                # Skip built-ins
                continue

            if name in self.__dict__:
                raise AttributeError(
                    f'Module at {path} cannot contain attributes {name} as it '
                    'overwrites an attribute of the same name in utils.params'
                )
            self.__setattr__(name, value)

params = Params()


def import_from_file(path):
    '''
    Return module object from a filepath
    '''
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} doesn\'t exist')
    
    name = path.split('/')[-1].split('.')[0]
    spec = spec_from_file_location(name, path)
    if spec is None:
        raise ValueError('Could not load module from "%s"' % path)
    
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m