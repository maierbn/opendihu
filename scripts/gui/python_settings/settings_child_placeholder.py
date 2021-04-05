from python_settings.settings_comment import SettingsComment

# a placeholder in a SettingsDict, which can be replaced with some SettingsDictEntrys
# this gets created when parsing python_options from possible_solver_combinations
# it is marked with the floating comment '### CHILD ###'
class SettingsChildPlaceholder(SettingsComment):
    def __init__(self, childnumber):
        super().__init__()
        self.comment = '### CHILD ' + str(childnumber) + ' ###'
        self.childnumber = int(childnumber)
