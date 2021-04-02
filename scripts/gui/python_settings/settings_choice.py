# holds 2 lists, one with default SettingsDictEntrys and one with alternative SettingsDictEntrys
class SettingsChoice:
    def __init__(self, defaults, alternatives):
        self.defaults = defaults
        self.alternatives = alternatives
