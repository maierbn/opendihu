import re

# import all settings-modules here, so we can only import this module to get them all
from python_settings.settings_activatable import *
from python_settings.settings_child_placeholder import *
from python_settings.settings_choice import *
from python_settings.settings_comment import *
from python_settings.settings_conditional import *
from python_settings.settings_container import *
from python_settings.settings_dict_entry import *
from python_settings.settings_empty_line import *
from python_settings.settings_list_entry import *


# this holds a complete settings.py by parsing its config-dict and storing the rest of the file in prefix and postfix
class PythonSettings():
    prefix = ''
    config_dict = None
    postfix = ''
    # takes a string of a settings.py and parses it

    def __init__(self, settings=None):
        if settings:
            # isolate content of config{} to settings and save the rest of the file settings_prefix and settings_postfix
            split1 = settings.split('config = {')
            self.prefix = split1[0][:-1]
            settings = split1[1]
            split2 = re.compile(r'(?m)^}').split(settings, 1)
            settings = split2[0]
            settings = '{' + settings + '}'
            self.postfix = split2[1][1:]

            # iterate over tokens to create SettingsDict
            self.config_dict = SettingsDict(settings)
            return None

    def __repr__(self):
        return self.prefix + '\nconfig = ' + str(self.config_dict) + self.postfix
