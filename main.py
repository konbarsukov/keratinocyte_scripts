from config.config_helper import Config
import ast

config = Config('config.cfg')

# print(config.project_name)
# print(config.global_name_list)
print(config.global_name_list('young_brd4_list'))
