import ConfigParser
import ast
import sys

class Config:

    def __init__(self, path='./config.cfg'):
        self.config = self.get_config(path)

        folders_dict = self.__create_section_dict('PROJECT_FOLDERS')
        for key, value in folders_dict.items():
            setattr(self, key, value)

    @property
    def project_name(self):
        return self.config.get('DEFAULT', 'project_name')

    @property
    def project_folder(self):
        return self.config.get('DEFAULT', 'project_folder')

    @property
    def pipeline_dir(self):
        return self.config.get('DEPENDENCIES', 'pipeline_dir')

    @property
    def genome(self):
        return self.config.get('DEPENDENCIES', 'genome')

    @property
    def annotation_file(self):
        return self.config.get('DEPENDENCIES', 'annotation_file')

    @property
    def mask_file(self):
        return self.config.get('DEPENDENCIES', 'mask_file')

    @property
    def genome_dir(self):
        return self.config.get('DEPENDENCIES', 'genome_directory')

    @property
    def gtf_file(self):
        return self.config.get('DEPENDENCIES', 'gtf_file')

    @property
    def py27_path(self):
        return self.config.get('DEPENDENCIES', 'py27_path')

    @property
    def crc_path(self):
        return self.config.get('DEPENDENCIES', 'crc_path')

    @property
    def fpkm_table(self):
        return self.config.get('DEPENDENCIES', 'fpkm_table')

    @property
    def data_tables_dict(self):
        return self.__create_section_dict('DATA_TABLES')

    # @property
    def global_name_list(self, key):

        v = self.__create_section_dict('GLOBAL_NAMES_LIST').get(key)
        if not v:
            return None

        try:
            list = ast.literal_eval(v)
        except SyntaxError:
            print('Error: {} can not be formated to list'.format(v))
            sys.exit(1)



        return ast.literal_eval(v)

    # __make_section_dict creates section dictionary
    # without DEFAULT options 'ProjectName' and 'ProjectFolder'
    def __create_section_dict(self, section_name):
        options_list = self.config.items(section_name)

        # make dictionary without default options
        default = ['project_name', 'project_folder']
        names_dict = {item[0]: item[1] for item in options_list if not default.count(item[0])}

        return names_dict

    @staticmethod
    def get_config(path):
        config = ConfigParser.ConfigParser()
        config.read(path)
        return config


if __name__ == "__main__":
    pass
