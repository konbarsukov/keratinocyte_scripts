import ConfigParser
import ast
import sys
import posixpath

'''
Config class contains all dependencies and methods to access them
'''
class Config:

    def __init__(self, path='./config.cfg'):
        self.config = self.get_config(path)

        folders_dict = self.__create_section_dict('PROJECT_FOLDERS')
        for key, value in folders_dict.items():
            setattr(self, key, value)

    '''
    Project name
    '''
    @property
    def project_name(self):
        return self.config.get('DEFAULT', 'project_name')

    '''
    Path to project folder
    '''
    @property
    def project_folder(self):
        path = self.config.get('DEFAULT', 'project_folder')
        _check_path(path)
        return path

    '''
    Path to pipeline_dir
    '''
    @property
    def pipeline_folder(self):
        path = self.config.get('DEPENDENCIES', 'pipeline_folder')
        _check_path(path)
        return path

    '''
    Genome name
    '''
    @property
    def genome_name(self):
        return self.config.get('DEPENDENCIES', 'genome_name')

    '''
    Path to annotation file
    '''
    @property
    def annotation_file(self):
        path = self.config.get('DEPENDENCIES', 'annotation_file')
        _check_path(path)
        return path

    '''
    Path to mask file
    '''
    @property
    def mask_file(self):
        path = self.config.get('DEPENDENCIES', 'mask_file')
        _check_path(path)
        return path

    '''
    Path to genome directory
    '''
    @property
    def genome_folder(self):
        path = self.config.get('DEPENDENCIES', 'genome_folder')
        _check_path(path)
        return path
    '''
    Path to gtf file
    '''
    @property
    def gtf_file(self):
        path = self.config.get('DEPENDENCIES', 'gtf_file')
        _check_path(path)
        return path

    '''
    Path to current python
    '''
    @property
    def py27_path(self):
        path = self.config.get('DEPENDENCIES', 'py27_path')
        _check_path(path)
        return path

    '''
    Path to crc binary file
    '''
    @property
    def crc_path(self):
        path = self.config.get('DEPENDENCIES', 'crc_path')
        _check_path(path)
        return path

    '''
    Path to dpkm_table (creates by cufflinks in step 0)
    '''
    @property
    def fpkm_table_file(self):
        path = self.config.get('DEPENDENCIES', 'fpkm_table_file')
        _check_path(path)
        return path

    '''
    Path to string_interaction_path
    '''
    @property
    def string_interaction_file(self):
        path = self.config.get('DEPENDENCIES', 'string_interaction_file')
        _check_path(path)
        return path

    '''
    Path to string_clustering_path
    '''
    @property
    def string_clustering_file(self):
        path = self.config.get('DEPENDENCIES', 'string_clustering_file')
        _check_path(path)
        return path

    '''
    get_data_table(self, key) returns data_tables or interrupts program with exit
    status 1
    '''
    def get_data_table(self, key):
        path = self.__get_dict('DATA_TABLES', key)
        _check_path(path)
        return path

    '''
    get_global_name(self, key) returns global_names or interrupts program with exit
    status 1
    '''
    def get_global_name(self, key):
        result = self.__get_dict('GLOBAL_NAMES_LIST', key)

        return result

    '''
    __make_section_dict() creates section dictionary
    without DEFAULT options 'ProjectName' and 'ProjectFolder
    '''
    def __create_section_dict(self, section_name):
        options_list = self.config.items(section_name)

        # make dictionary without default options
        default = ['project_name', 'project_folder']
        names_dict = {item[0]: item[1] for item in options_list if not default.count(item[0])}

        return names_dict

    '''
    __get_dict() returns section dictionary
    if key not found method interrupts program with status 1
    if key value can not be converted to python literal,
    method interrupts program with status 1
    '''
    def __get_dict(self, section, key):
        result = None
        v = self.__create_section_dict(section).get(key)
        if not v:
            print('Config: in section {} value {} for key {} not found'.format(section, v, key))
            sys.exit(1)
        try:
            result = ast.literal_eval(v)
        except SyntaxError:
            print('Config: in section {} value {} for key {} can not be formated'.format(section, v, key))
            sys.exit(1)

        return result

    '''
    get_config() return initialized config object
    '''
    @staticmethod
    def get_config(path):
        _check_path(path)
        config = ConfigParser.ConfigParser()
        config.read(path)
        return config

'''
_check_path() closes program if file is not existed or exit with status 1
'''
def _check_path(path):
    if not posixpath.exists(path):
        print('Config: file or dir {} is not existed'.format(path))
        sys.exit(1)
    print('Config: found file or dir {}'.format(path))



if __name__ == "__main__":
    pass
