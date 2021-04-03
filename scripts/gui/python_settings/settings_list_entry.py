from python_settings.settings_activatable import Activatable

# normal entry for a SettingsList
class SettingsListEntry(Activatable):
    def __init__(self, value=None, comment=None):
        self.value = value
        self.comments = []
        if comment:
            self.comments.append('#' + comment)
