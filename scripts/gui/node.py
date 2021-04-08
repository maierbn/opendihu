import copy

from helpers import printe, indent, Error, Info, Warning
from python_settings.python_settings import *

class Childs():
    def __init__(self, node):
        self.node = node
        # node.name is not populated yet, so we can't populate self.__childs with placeholders
        self.populated = False

    def populate(self):
        self.populated = True

        self.__childs = []
        self.combinations = self.node.combinations

        if not isinstance(self.node, PlaceholderNode) and self.node.name in self.node.combinations and "template_arguments" in self.node.combinations[self.node.name]:
            template_arguments = self.node.combinations[self.node.name]["template_arguments"]
            childs_count = len(template_arguments)
            self.__childs = [None] * childs_count
            if "template_arguments_needed" in self.node.combinations[self.node.name]:
                self.childs_count_needed = self.node.combinations[
                    self.node.name]["template_arguments_needed"]
            else:
                self.childs_count_needed = childs_count
            for i in range(childs_count):
                self.__add_placeholder_i(i)

    def __add_placeholder_i(self, i):
        if i < self.childs_count_needed:
            self.__childs[i] = (PlaceholderNode(
                self.combinations, needed=True))
        else:
            self.__childs[i] = (PlaceholderNode(
                self.combinations, needed=False))
        self.__childs[i].parent = self.node

    def get_childs(self):
        if not self.populated:
            self.populate()
        return self.__childs

    def get_real_childs(self):
        if not self.populated:
            self.populate()
        ret = []
        for child in self.__childs:
            if not isinstance(child, PlaceholderNode):
                ret.append(child)
        return ret

    def delete(self, child):
        if not self.populated:
            self.populate()
        for i in range(len(self.__childs)):
            if self.__childs[i] == child:
                self.__add_placeholder_i(i)

    def replace(self, child_old, child_new):
        if not self.populated:
            self.populate()
        child_new.parent = self.node
        for i in range(len(self.__childs)):
            if self.__childs[i] == child_old:
                self.__childs[i] = child_new
                return
        #printe('failed to replace node')

    # normally gets called while parsing cpp-code
    def replace_next_placeholder(self, child):
        if not self.populated:
            self.populate()
        for c in self.get_childs():
            if isinstance(c, PlaceholderNode):
                self.replace(c, child)
                return
        # force adding the child if no PlaceholderNodes are left in self.__childs (in case of unknown templates)
        child.parent = self.node
        self.__childs.append(child)

    def clear(self):
        self.populate()

# this class represents a Node in the structure tree (Example.root e.g. is such a Node)
class Node:
    def __init__(self, combinations):
        self.combinations = combinations
        self.name = ''
        self.comment = ''
        self.can_have_childs = False
        self.childs = Childs(self)

        self.parent = None

        self.settings_dict = SettingsDict()
        self.settings_container_default = None

    def get_int_replacement(self, integer):
        replacement = Node(self.combinations)
        replacement.name = integer
        return replacement

    def get_contextual_description(self):
        # special case for rootnode
        if not self.parent:
            return ''
        try:
            for i in range(len(self.parent.childs.get_childs())):
                if self == self.parent.childs.get_childs()[i]:
                    child_index = i
                    break
            (node_description,
             _possible_node_names) = self.combinations[self.parent.name]["template_arguments"][child_index]
            return node_description
        except:
            return 'UNKNOWN'

    def get_possible_replacements(self):
        # special case for rootnode
        if not self.parent:
            return []
        for i in range(len(self.parent.childs.get_childs())):
            if self == self.parent.childs.get_childs()[i]:
                child_index = i
                break
        (possible_node_description,
         possible_node_names) = self.combinations[self.parent.name]["template_arguments"][child_index]
        possible_replacements = []
        for name in possible_node_names:
            possible_replacement = Node(self.combinations)
            possible_replacement.name = name
            # possible_replacement.add_missing_default_python_settings()
            try:  # try because of Integer name not found in self.combinations
                if "template_arguments" in self.combinations[name]:
                    possible_replacement.can_have_childs = True
            except:
                pass
            possible_replacements.append(possible_replacement)
        # TODO sort by occurence in examples
        return (possible_node_description, possible_replacements)

    # this is not in __init__(), because self.name (used here) gets defined later
    def get_default_python_settings_dict(self):
        if self.settings_container_default == None:
            try:
                self.settings_container_default = self.combinations[self.name]["python_options"]
            except:
                #printe('no python_options found for ' + str(self.name))
                # return None if nothing found
                self.settings_container_default = None
        return self.settings_container_default

    # returns self.settings_dict with SettingsChildPlaceholders replaced with child dicts
    def get_python_settings_dict_recursive(self):
        # deepcopy self.settings_dict so we don't replace the SettingsChildPlaceholders in it
        own_dict = copy.deepcopy(self.settings_dict)
        for child in self.childs.get_real_childs():
            #print('child: ' + child.name)
            # for every child replace the ### CHILD ### placeholder with the childs dict
            try:
                child_dict = child.get_python_settings_dict_recursive()
                #print('replacing ' + child.name + ' placeholder with:\n' + str(child_dict))
                own_dict.replaceChildPlaceholder(child_dict)
            except:
                pass
                #print('failed to replace SettingsChildPlaceholder ' + child.name)
        return own_dict

    # inserts deactivated placeholders for all possible settings into self.settings_dict
    def insert_missing_default_python_settings_deactivated(self, settings_global_dict, recurse_childs=True, self_settings_container=None, settings_container_default=None):
        # counter for inserted settings
        changes = 0
        recurse_childs_at_the_end = False

        # init stuff and recurse childs if we are on the outer level of a node
        if self_settings_container == None and settings_container_default == None:
            # init settings_container_default
            settings_container_default = self.get_default_python_settings_dict()
            if settings_container_default == None:
                # TODO maybe add error (no defaults found for self.name)
                return 0

            # init self_settings_container
            if self.settings_dict == None:
                self.settings_dict = SettingsDict()
            self_settings_container = self.settings_dict

            self.childs_with_placeholders = []
            # recurse childs
            if recurse_childs:
                recurse_childs_at_the_end = True

        # handle SettingsList
        if isinstance(self_settings_container, SettingsList) and isinstance(settings_container_default, SettingsList):
            # if list is empty -> add all list entries + list_comprehension
            if len(self_settings_container) == 0 and len(settings_container_default) > 0:
                self_settings_container.list_comprehension = settings_container_default.list_comprehension
                for settings_container_default_entry in settings_container_default:
                    default_entry = copy.deepcopy(settings_container_default_entry)
                    if not isinstance(default_entry.value, str):
                        default_entry.value = type(default_entry.value)()
                    # deactivate added entries
                    default_entry.activated = False
                    self_settings_container.append(default_entry)
            ## recurse
            # if we have multiple entries just use settings_container_default[0] for all of them
            for entry in self_settings_container:
                # add parent_info
                entry.parent = self_settings_container

                if isinstance(entry, SettingsListEntry):
                    if not isinstance(entry.value, str):
                        # add parent_info
                        entry.value.parent = entry

                        settings_container_default_recurse = None
                        try:
                            settings_container_default_recurse = settings_container_default[0].value
                        except: pass
                        changes = changes + self.insert_missing_default_python_settings_deactivated(self_settings_container=entry.value, recurse_childs=recurse_childs,
                                                                                     settings_container_default=settings_container_default_recurse, settings_global_dict=settings_global_dict,)
        # handle SettingsDict
        elif isinstance(self_settings_container, SettingsDict) and isinstance(settings_container_default, SettingsDict):
            # resolve all SettingsChoice
            if isinstance(settings_container_default, SettingsDict):
                settings_container_default_resolved = SettingsDict()
                for entry in settings_container_default:
                    if isinstance(entry, SettingsChoice):
                        settings_container_default_resolved.extend(entry.defaults)
                        settings_container_default_resolved.extend(entry.alternatives)
                    else:
                        settings_container_default_resolved.append(entry)
                settings_container_default = settings_container_default_resolved

            # add missing default settings to this level (NOT recursive)
            for entry in settings_container_default:
                if isinstance(entry, SettingsChildPlaceholder):
                    self.childs_with_placeholders.append(self.childs.get_childs()[entry.childnumber])
                    placeholder_already_added = False
                    for e in self_settings_container:
                        if isinstance(e, SettingsChildPlaceholder) and e.childnumber == entry.childnumber:
                            placeholder_already_added = True
                            break
                    if not placeholder_already_added:
                        self_settings_container.append(entry)

                elif isinstance(entry, SettingsDictEntry):
                    # add default-entry if we don't have the key already
                    if not self_settings_container.has_key(entry.key):
                        entry_copy = copy.deepcopy(entry)
                        if not isinstance(entry_copy.value, str):
                            # if entry_copy.value is not a string (then it is SettingsDict or SettingsList) -> replace it with empty one
                            entry_copy.value = type(entry_copy.value)()
                        else:
                            changes = changes + 1
                            #print('added: ' + entry_copy.key + ':' + entry_copy.value)
                        # deactivate entry before adding
                        entry_copy.activated = False
                        self_settings_container.append(entry_copy)

                elif isinstance(entry, SettingsMesh) or isinstance(entry, SettingsSolver):
                    #local
                    dict_to_append_to = self_settings_container
                    # add all missing keys recursively
                    changes = changes + self.insert_missing_default_python_settings_deactivated(
                        self_settings_container=dict_to_append_to, recurse_childs=recurse_childs, settings_container_default=entry, settings_global_dict=settings_global_dict)
                    # global
                    # add e.g. Meshes : {} if not there
                    if not settings_global_dict.has_key(entry.global_key):
                        settings_global_dict.append(
                            SettingsDictEntry(entry.global_key, SettingsDict()))
                    dict = settings_global_dict.get_value(entry.global_key)
                    if not isinstance(dict, SettingsDict):
                        #printe('we have to add to global, but global is not a dict (propably it is a variable, we cannot add to)')
                        continue
                    # get the name e.g. mesh0
                    if self_settings_container.has_key(entry.name_key):
                        name = self_settings_container.get_value(
                            entry.name_key)
                    else:
                        i = 0
                        name = ''
                        while True:
                            name = '"' + entry.name_prefix + str(i) + '"'
                            if not dict.has_key(name):
                                break
                            i = i+1
                        name_entry = SettingsDictEntry(entry.name_key, name)
                        name_entry.activated = False
                        self_settings_container.append(name_entry)
                        #changes = changes + 1
                    # add global dict entry if not there (e.g. Meshes : { mesh0 : {} })
                    if not dict.has_key(name):
                        global_dict_entry = SettingsDictEntry(name, SettingsDict())
                        global_dict_entry.activated = False
                        dict.append(global_dict_entry)
                        #changes = changes + 1
                    dict_to_append_to = dict.get_value(name)

                    # add all missing keys recursively
                    changes = changes + self.insert_missing_default_python_settings_deactivated(
                        self_settings_container=dict_to_append_to, recurse_childs=recurse_childs, settings_container_default=entry, settings_global_dict=settings_global_dict)

            # recurse levels
            for entry in self_settings_container:
                # add parent_info
                entry.parent = self_settings_container

                if isinstance(entry, SettingsDictEntry):
                    # add doc_link to recursive SettingsDictEntrys
                    # and add default_comment (used in change-settings-gui)
                    try:
                        default_entry = settings_container_default.get_entry(entry.key)
                        # the order is important here, because accessing comments[0] could fail
                        entry.doc_link = default_entry.doc_link
                        entry.default_comment = default_entry.comments[0]
                    except: pass

                    # recurse all keys that are no strings, those are SettingsDict and SettingsList
                    if not isinstance(entry.value, str) and settings_container_default.has_key(entry.key):
                        # add parent_info
                        entry.value.parent = entry

                        settings_container_default_recurse = settings_container_default.get_value(
                            entry.key)
                        changes = changes + self.insert_missing_default_python_settings_deactivated(
                            self_settings_container=entry.value, recurse_childs=recurse_childs, settings_container_default=settings_container_default_recurse, settings_global_dict=settings_global_dict)

        if recurse_childs_at_the_end:
            # only recurse childs, for that we have seen or set a child_placeholders
            # if e.g. a dict in which a SettingsChildPlaceholder should be is set to an external variable,
            # we do not want to add those settings again to the childnode
            childs_recurse = []
            for child in self.childs_with_placeholders:
                # ignore childs that have an Integer as name
                try:
                    int(child.name)
                except:
                    childs_recurse.append(child)
            for child in childs_recurse:
                changes = changes + child.insert_missing_default_python_settings_deactivated(
                    settings_global_dict=settings_global_dict, recurse_childs=recurse_childs)

        # except:
        #    printe('something went wrong while adding missing python-settings')

        return changes


    # activates all deactivated default python-settings, respecting SettingsChoice (only activates SettingsChoice.defaults)
    def activate_all_default_python_settings(self, settings_global_dict, recurse_childs=True, self_settings_container=None, settings_container_default=None):
        # counter for added settings
        changes = 0
        recurse_childs_at_the_end = False

        # init stuff and recurse childs if we are on the outer level of a node
        if self_settings_container == None and settings_container_default == None:
            # init settings_container_default
            settings_container_default = self.get_default_python_settings_dict()
            if settings_container_default == None:
                return 0

            # init self_settings_container
            if self.settings_dict == None:
                self.settings_dict = SettingsDict()
            self_settings_container = self.settings_dict

            if recurse_childs:
                recurse_childs_at_the_end = True

        # handle SettingsList
        if isinstance(self_settings_container, SettingsList) and isinstance(settings_container_default, SettingsList):
            ## recurse
            # if we have multiple entries just use settings_container_default[0] for all of them
            for entry in self_settings_container:
                if isinstance(entry, SettingsListEntry):
                    if not isinstance(entry.value, str):
                        settings_container_default_recurse = None
                        try:
                            settings_container_default_recurse = settings_container_default[0].value
                        except: pass
                        changes = changes + self.activate_all_default_python_settings(self_settings_container=entry.value, recurse_childs=recurse_childs,
                                                                                     settings_container_default=settings_container_default_recurse, settings_global_dict=settings_global_dict,)
        # handle SettingsDict
        elif isinstance(self_settings_container, SettingsDict) and isinstance(settings_container_default, SettingsDict):
            # resolve all SettingsChoice to defaults
            if isinstance(settings_container_default, SettingsDict):
                settings_container_default_resolved = SettingsDict()
                for entry in settings_container_default:
                    if isinstance(entry, SettingsChoice):
                        add_alternatives = False
                        # look if some settings from alternatives already exist AND are activated, in that case add alternatives instead of defaults
                        for e in entry.alternatives:
                            activated = True
                            try:
                                activated = self_settings_container.get_entry(e.key).activated
                            except: pass
                            if self_settings_container.has_key(e.key) and activated:
                                add_alternatives = True
                                break
                        if add_alternatives:
                            for e in entry.alternatives:
                                settings_container_default_resolved.append(e)
                        else:
                            for e in entry.defaults:
                                settings_container_default_resolved.append(e)
                    else:
                        settings_container_default_resolved.append(entry)
                settings_container_default = settings_container_default_resolved

                for entry in settings_container_default:
                    if isinstance(entry, SettingsDictEntry):
                        # activate entries
                        # try because e.g. the key "meshName" is not found in self_settings_container
                        try:
                            e = self_settings_container.get_entry(entry.key)
                            if e.activated == False:
                                e.activated = True
                                changes = changes + 1
                        except: pass

                    elif isinstance(entry, SettingsMesh) or isinstance(entry, SettingsSolver):
                        dict_to_append_to = self_settings_container
                        name_entry = self_settings_container.get_entry(entry.name_key)
                        if name_entry.activated:
                            activate_global = True
                        else:
                            # resolve SettingsChoice inside Mesh/Solver
                            entry_resolved = SettingsDict()
                            for e in entry:
                                if isinstance(e, SettingsChoice):
                                    entry_resolved.extend(e.defaults)
                                    entry_resolved.extend(e.alternatives)
                                else:
                                    entry_resolved.append(e)

                            # check if any of the entries is activated, if not -> global by default
                            activate_global = True
                            for e in entry_resolved:
                                if self_settings_container.get_entry(e.key).activated:
                                    activate_global = False
                                    break

                        if activate_global:
                            # activate global
                            name_entry.activated = True
                            meshes = settings_global_dict.get_value(entry.global_key)
                            if not isinstance(meshes, SettingsDict):
                                printe('we have to activate keys in global, but global is not a dict (propably it is a variable, we cannot add to)')
                                continue
                            global_dict_entry = meshes.get_entry(name_entry.value)
                            global_dict_entry.activated = True
                            dict_to_append_to = global_dict_entry.value
                        else:
                            # activate local
                            dict_to_append_to = self_settings_container
                        # activate keys recursively
                        changes = changes + self.activate_all_default_python_settings(
                            self_settings_container=dict_to_append_to, recurse_childs=recurse_childs, settings_container_default=entry, settings_global_dict=settings_global_dict)

            for entry in self_settings_container:
                # recurse levels
                if isinstance(entry, SettingsDictEntry):
                    # recurse all keys that are no strings, those are SettingsDict and SettingsList
                    if not isinstance(entry.value, str) and settings_container_default.has_key(entry.key):
                        settings_container_default_recurse = settings_container_default.get_value(
                            entry.key)
                        changes = changes + self.activate_all_default_python_settings(
                            self_settings_container=entry.value, recurse_childs=recurse_childs, settings_container_default=settings_container_default_recurse, settings_global_dict=settings_global_dict)

        if recurse_childs_at_the_end:
            childs_recurse = []
            for child in self.childs.get_real_childs():
                # ignore childs that have an Integer as name
                try:
                    int(child.name)
                except:
                    childs_recurse.append(child)
            for child in childs_recurse:
                changes = changes + child.activate_all_default_python_settings(
                    settings_global_dict=settings_global_dict, recurse_childs=recurse_childs)

        return changes


    # delete self.settings_dict recursively
    def delete_python_settings_recursive(self):
        self.settings_dict = SettingsDict()
        for child in self.childs.get_real_childs():
            child.delete_python_settings_recursive()

    # parse PythonSettings and keep prefix and postfix
    # returns a list of Warnings
    def parse_python_settings(self, python_settings, keep_entries_that_have_no_default, recurse_childs=True):
        # remove all old python-settings
        if recurse_childs:
            self.delete_python_settings_recursive()
        else:
            self.settings_dict = SettingsDict()
        self.settings_dict_prefix = python_settings.prefix
        self.settings_dict_postfix = python_settings.postfix
        (_, warnings) = self.parse_python_settings_recursive(python_settings.config_dict,
                                                             keep_entries_that_have_no_default=keep_entries_that_have_no_default, recurse_childs=recurse_childs, warnings=[])
        warnings.append(Info('written python-settings for ' + str(self.name)))
        return warnings

    # parse a python_settings_dict and add it to this Node and its childs
    # returns (rest, warnings)
    # rest is only not None in recursive calls on childs and can be ignored if called from outside
    # warnings is a list of Warnings
    def parse_python_settings_recursive(self, settings_container, keep_entries_that_have_no_default=False, self_settings_container=None, settings_container_default=None, is_called_on_child=False, warnings=[], recurse_childs=True):
        # TODO SettingsConditional
        # TODO if a SettingsDictEntry has no comment, add the default-comment (extra function)
        # these should be given as parameters when recursion happens on the same object again
        # self_settings_container is a reference to the sub-SettingsDict in self.settings_dict we are currently handling
        # settings_container_default is a reference to the sub-SettingsDict in settings_container_default we are currently handling
        if self_settings_container == None:
            if self.settings_dict == None:
                self.settings_dict = SettingsDict()
            self_settings_container = self.settings_dict
        if settings_container_default == None:
             settings_container_default = self.get_default_python_settings_dict()

        # create a list with all child_placeholders at this level
        child_placeholders = []
        # here we assume that SettingsChildPlaceholders never accur inside SettingsChoice (SettingsChoice of settings_container_default are not resolved here)
        for entry in settings_container_default:
            if isinstance(entry, SettingsChildPlaceholder):
                child_placeholders.append(entry)

        # if is_called_on_child:
        rest = SettingsDict()

        for i in range(len(settings_container)):
            entry = settings_container[i]
            #try:
            #    print(type(entry.value))
            #    print(entry.value)
            #except: pass

            # always keep Solvers and Meshes and meta
            if isinstance(entry, SettingsDictEntry) and (entry.key == '"Solvers"' or entry.key == '"Meshes"' or entry.key == '"meta"' or entry.key == '"MappingsBetweenMeshes"'):
                self_settings_container.append(entry)
                continue

            if isinstance(entry, SettingsDictEntry) or isinstance(entry, SettingsListEntry):
                # resolve SettingsMesh and SettingsSolver
                # we do this before resolving SettingsChoice to enable resolving SettingsChoice inside SettingsMesh and SettingsSolver
                # TODO maybe its needed to resolve choices before and after resolving SettingsMesh and SettingsSolver
                # TODO either accept local or global defs
                if isinstance(settings_container_default, SettingsDict):
                    settings_container_default_resolved = SettingsDict()
                    for e in settings_container_default:
                        if isinstance(e, SettingsMesh) or isinstance(e, SettingsSolver):
                            #print(e.name_key)
                            settings_container_default_resolved.append(
                                SettingsDictEntry(e.name_key, ""))
                            for default_e in e:
                                settings_container_default_resolved.append(
                                    default_e)
                        else:
                            settings_container_default_resolved.append(e)
                    settings_container_default = settings_container_default_resolved

                # resolve SettingsChoice inside settings_container_default
                # TODO maybe do this a little bit smarter (do not accept all things from the SettingsChoice (either defaults or alternatives))
                if isinstance(settings_container_default, SettingsDict):
                    settings_container_default_resolved = SettingsDict()
                    for e in settings_container_default:
                        if isinstance(e, SettingsChoice):
                            for default_e in e.defaults:
                                settings_container_default_resolved.append(
                                    default_e)
                            for alternative_e in e.alternatives:
                                settings_container_default_resolved.append(
                                    alternative_e)
                        else:
                            settings_container_default_resolved.append(e)
                    settings_container_default = settings_container_default_resolved

                # check if the entry is a SettingsListEntry or if it is a SettingsDictEntry and its key is in the defaults
                if isinstance(entry, SettingsListEntry) or (isinstance(settings_container_default, SettingsDict) and settings_container_default.has_key(entry.key)):
                    # key exists in defaults
                    # if the value is a SettingsContainer -> recurse
                    if isinstance(entry.value, SettingsContainer):
                        # create a new entry
                        settings_container_default_recurse = []
                        if isinstance(entry, SettingsDictEntry):
                            new_entry = SettingsDictEntry()
                            new_entry.key = entry.key
                            settings_container_default_recurse = settings_container_default.get_value(
                                entry.key)
                        else:
                            new_entry = SettingsListEntry()
                            # TODO here we assume that there is only one SettingsListEntry in python_options we just use the first one
                            try:
                                settings_container_default_recurse = settings_container_default.get_first_SettingsListEntry().value
                            except: pass
                        new_entry.comments = entry.comments
                        ## add doc_link to recursive SettingsDictEntrys
                        ## and add default_comment (used in change-settings-gui)
                        #try:
                        #    default_entry = settings_container_default.get_entry(entry.key)
                        #    # the order is important here, because accessing comments[0] could fail
                        #    new_entry.doc_link = default_entry.doc_link
                        #    new_entry.default_comment = default_entry.comments[0]
                        #except: pass

                        # add an empty list or dict to the new entry
                        if isinstance(entry.value, SettingsDict):
                            new_entry.value = SettingsDict()
                        else:
                            new_entry.value = SettingsList()
                            new_entry.value.list_comprehension = entry.value.list_comprehension
                        # add the new entry
                        self_settings_container.append(new_entry)

                        # recurse the new entry
                        (_, warnings) = self.parse_python_settings_recursive(entry.value, self_settings_container=new_entry.value, settings_container_default=settings_container_default_recurse,
                                                                             keep_entries_that_have_no_default=keep_entries_that_have_no_default, warnings=warnings, recurse_childs=recurse_childs)
                    else:
                        # if the entry.value is no SettingsContainer -> just append the entry
                        ## add doc_link
                        #try:
                        #    default_entry = settings_container_default.get_entry(entry.key)
                        #    entry.doc_link = default_entry.doc_link
                        #except: pass
                        self_settings_container.append(entry)
                else:
                    # entry is a SettingsDictEntry
                    # and the key does not exist in defaults
                    # this entry possibly belongs to a child
                    new_dict = SettingsDict()
                    new_dict.append(entry)
                    # try giving the entry to the childs
                    if recurse_childs:
                        for j in range(len(child_placeholders)):
                            #print('trying to give ' + str(entry.key) + ' to child ' + str(child_placeholders[j].childnumber))
                            child = self.childs.get_childs(
                            )[child_placeholders[j].childnumber]
                            if not isinstance(child, PlaceholderNode):
                                (new_dict, warnings) = child.parse_python_settings_recursive(
                                    new_dict, keep_entries_that_have_no_default=keep_entries_that_have_no_default, is_called_on_child=True, recurse_childs=recurse_childs, warnings=warnings)
                            if len(new_dict) == 0:
                                break

                    if len(new_dict) > 0:
                        # this entry is not in defaults and does not belong to a child
                        if is_called_on_child:
                            rest.append(entry)
                        elif keep_entries_that_have_no_default:
                            entry.is_unknown = True
                            warnings.append(
                                Warning(entry.key + ' is an unknown setting -> added it to ' + str(self.name) + ' anyway'))
                            self_settings_container.append(entry)
                        else:
                            warnings.append(
                                Warning(entry.key + ' is an unknown setting -> it was NOT added to ' + str(self.name)))

            elif isinstance(entry, SettingsChildPlaceholder):
                # don't add placeholders found in SettingsDict
                pass
            else:
                # always add SettingsComments
                self_settings_container.append(entry)

        # insert all SettingsChildPlaceholders that are in python_options on this level
        # we do this here, to always have the childs on the bottom for consistency
        # TODO maybe append them to self_settings_container after first use and only the not used ones here at the bottom
        for child_placeholder in child_placeholders:
            self_settings_container.append(child_placeholder)

        return (rest, warnings)

    # this function converts the tree under this Node to a pretty string
    # you can print the string to visualize the tree
    # this is also used to created cpp-source-code from a tree
    def __repr__(self):
        comment = ''
        if self.comment != '':
            comment = '//' + self.comment + '\n'
        return comment + self.repr_recursive(0)

    def repr_recursive(self, depth):
        indentation = '  ' * depth
        indentation_child = '  ' * (depth + 1)
        childs_string = ''
        for i in range(len(self.childs.get_real_childs())):
            child = self.childs.get_real_childs()[i]
            #print(child.name + ' : ' + child.comment)
            comment_string = ''
            if child.comment != '':
                comment_string = ' //' + child.comment
            # add ',' if this is not the last child
            if i < len(self.childs.get_real_childs()) - 1:
                comment_string = ',' + comment_string
            if childs_string == '':
                childs_string = '\n' + indentation_child + \
                    child.repr_recursive(depth + 1) + comment_string
            else:
                childs_string = childs_string + '\n' + indentation_child + \
                    child.repr_recursive(depth + 1) + comment_string
        if childs_string == '':
            if self.can_have_childs:
                return self.name + '<>'
            else:
                return self.name
        else:
            return self.name + '<' + childs_string + '\n' + indentation + '>'

    # this function can compare the cpp of this Node to another Node
    # returning True if equal, False otherwise
    def compare_cpp(self, node):
        if self.name != node.name:
            return False
        if self.comment != node.comment:
            #print(self.name + ' : ' + self.comment + ' =! ' + node.comment)
            return False
        if self.can_have_childs != node.can_have_childs:
            return False
        if len(self.childs.get_real_childs()) != len(node.childs.get_real_childs()):
            return False
        for i in range(len(self.childs.get_real_childs())):
            if self.childs.get_real_childs()[i].compare_cpp(node.childs.get_real_childs()[i]) == False:
                return False
        return True

    # returns list of Errors or empty list if the tree is a valid combination of templates (non recursive)
    def validate_cpp_src(self):
        res = []
        # if not RootNode:
        if self.parent:
            try:
                p_wanted_childs = self.combinations[self.parent.name].get(
                    "template_arguments", [])
                p_argument_count_max = len(p_wanted_childs)
            except:
                res.append(Error(str(self.parent.name) + ' is unknown'))
                return res
            i = 0
            for child in self.parent.childs.get_childs():
                if child == self:
                    break
                i = i + 1

            (_template_argument_description,
             possible_template_arguments) = p_wanted_childs[i]
            if i > p_argument_count_max:
                res.append(Error(str(self.parent.name) + ' only accepts ' +
                                 str(p_argument_count_max) + ' template_arguments'))
            else:
                if possible_template_arguments == ["Integer"]:
                    try:
                        int(self.name)
                    except:
                        res.append(
                            Error(str(self.name) + ' is not an Integer'))
                elif self.name not in possible_template_arguments:
                    res.append(Error(str(self.name) + ' is not in the list of possible template_arguments for ' +
                                     self.parent.name + '\n' + 'possible template_arguments are: ' + str(possible_template_arguments)))

        return res

    # returns list of Errors or empty list if the tree is a valid combination of templates
    def validate_cpp_src_recursive(self):
        return self._validate_cpp_src_recursive(res=[])

    # helper function
    def _validate_cpp_src_recursive(self, recurse=True, node=None, res=[]):
        if not node:
            node = self
        try:
            wanted_childs = self.combinations[node.name].get(
                "template_arguments", [])
            argument_count_max = len(wanted_childs)
            argument_count_min = self.combinations[node.name].get(
                "template_arguments_needed", argument_count_max)
        except:
            # if the key node.name does not exist, we are at the bottom
            return res
        child_count = len(node.childs.get_real_childs())
        if child_count < argument_count_min:
            res.append(Error(str(node.name) + ' needs at least ' + str(argument_count_min) +
                             ' template_arguments but only ' + str(child_count) + ' template_arguments given'))
        if child_count > argument_count_max:
            res.append(Error(str(node.name) + ' only accepts ' + str(argument_count_max) +
                             ' template_arguments ' + str(child_count) + ' template_arguments given'))
        for i in range(len(node.childs.get_real_childs())):
            try:
                (_template_argument_description,
                 possible_template_arguments) = wanted_childs[i]
                if possible_template_arguments == ["Integer"]:
                    try:
                        int(node.childs.get_real_childs()[i].name)
                    except:
                        res.append(
                            Error(str(node.childs.get_real_childs()[i].name) + ' is not an Integer'))
                elif node.childs.get_real_childs()[i].name not in possible_template_arguments:
                    res.append(Error(str(node.childs.get_real_childs()[
                               i].name) + ' is not in the list of possible template_arguments for ' + node.name + '\n' + 'possible template_arguments are: ' + str(possible_template_arguments)))
            except:
                pass
            if recurse:
                res = self._validate_cpp_src_recursive(
                    node=node.childs.get_real_childs()[i], res=res)
        return res

class PlaceholderNode(Node):
    def __init__(self, combinations, needed):
        self.needed = needed
        super().__init__(combinations)
