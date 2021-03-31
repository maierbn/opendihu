from node import Node
from helpers import printe, indent, Error, Info
from python_settings.python_settings import *

# specialized Node with some extra stuff


class RootNode(Node):
    def __init__(self, combinations):
        super().__init__(combinations)

        self.name = 'GLOBAL'

        self.settings_dict_prefix = ''
        self.settings_dict_postfix = ''

        self.insert_missing_default_python_settings_deactivated(self.settings_dict)
        self.activate_all_default_python_settings(self.settings_dict)

    def get_python_settings(self):
        settings_dict = self.get_python_settings_dict_recursive()
        if not settings_dict:
            settings_dict = SettingsDict()
        python_settings = PythonSettings()
        python_settings.config_dict = settings_dict
        python_settings.prefix = self.settings_dict_prefix
        python_settings.postfix = self.settings_dict_postfix
        return python_settings
