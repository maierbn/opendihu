#!/usr/bin/env python3
import sys
import subprocess
import os
import copy

import gi
gi.require_version('Gtk', '3.0')
gi.require_version('GtkSource', '3.0')
from gi.repository import Gtk, Gio, GtkSource, GObject, Gdk, GLib

from cpp_tree import CPPTree
from root_node import RootNode
from node import PlaceholderNode
from python_settings.python_settings import *
from helpers import Message, Error, Info, Warning


# stores a node + its depth to view it in a ListBox
class NodeLine(GObject.GObject):
    def __init__(self, node, depth):
        GObject.GObject.__init__(self)
        self.node = node
        self.depth = depth

class SettingsLine(GObject.GObject):
    def __init__(self, settings, depth, parent):
        GObject.GObject.__init__(self)
        self.settings = settings
        self.depth = depth
        self.parent = parent


class ListBoxRowWithNode(Gtk.ListBoxRow):
    def add_node(self, node):
        self.node = node

class ListBoxRowWithSettings(Gtk.ListBoxRow):
    def add_node(self, settings):
        self.settings = settings


class DiscardNodeChangesDialog(Gtk.Dialog):
    def __init__(self, parent):
        Gtk.Dialog.__init__(self, title="Warning",
                            transient_for=parent, flags=0)
        self.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                         'discard changes', Gtk.ResponseType.OK)
        self.set_default_size(150, 100)
        label = Gtk.Label(
            label="Warning: Some changes to the current Node are not applied yet, and will be discarded!")
        box = self.get_content_area()
        box.add(label)
        self.show_all()

class PythonSettingsChangeWindow(Gtk.Window):
    def __init__(self, settings, main_window):
        super(PythonSettingsChangeWindow, self).__init__()

        grid = Gtk.Grid()
        self.add(grid)

        if isinstance(settings, SettingsDictEntry):
            title = settings.key
        elif isinstance(settings, SettingsListEntry):
            title = "list-entry"
        else:
            self.close()
        label_title = Gtk.Label(label=str(title))
        grid.add(label_title)

        language_manager = GtkSource.LanguageManager()
        text_view_python_code = GtkSource.View()
        text_view_python_code.get_buffer().set_language(
            language_manager.get_language('python3'))
        text_view_python_code.set_vexpand(True)
        text_view_python_code.set_hexpand(True)
        scroll_python_code = Gtk.ScrolledWindow()
        scroll_python_code.add(text_view_python_code)
        scroll_python_code.set_min_content_height(200)
        scroll_python_code.set_min_content_width(400)
        grid.attach_next_to(scroll_python_code, label_title,
                            Gtk.PositionType.BOTTOM, 1, 1)

        text = str(settings.value)
        text_view_python_code.get_buffer().set_text(text)
        text_original = copy.deepcopy(text)

        button_grid = Gtk.Grid()
        grid.attach_next_to(button_grid, scroll_python_code,
                            Gtk.PositionType.BOTTOM, 1, 1)

        def on_button_ok(_):
            text_bounds = text_view_python_code.get_buffer().get_bounds()
            text = text_view_python_code.get_buffer().get_text(text_bounds[0], text_bounds[1], True)
            if text != text_original:
                main_window.cpp_tree.undo_stack.duplicate_current_state()
                # activate settings recurse up
                settings.activate_recursive()

                # this is a little hacky:
                # we first set the value to the string
                settings.value = text
                # then we reparse the whole python_settings string
                new_settings_str = str(main_window.cpp_tree.get_python_settings())
                rets = main_window.cpp_tree.parse_python_settings(new_settings_str, undoable=False)
                main_window.log_append_message(rets)
                if not isinstance(rets, Error):
                    main_window.redraw_python()
                    self.close()
                else:
                    # revert to the undo point, because changes were syntactically incorrect
                    main_window.cpp_tree.undo_stack.undo()
                    main_window.cpp_tree.undo_stack.remove_future()
                    # we have to redraw_all after undo
                    main_window.redraw_all()
                    # if we don't close here settings is a dangling pointer to an old root, because of the undo
                    self.close()
            else:
                # nothing changed
                self.close()

        def on_button_cancel(_):
            self.close()

        button_ok = Gtk.Button(label='ok')
        button_ok.set_hexpand(True)
        button_ok.connect("clicked", on_button_ok)
        button_grid.add(button_ok)

        button_cancel = Gtk.Button(label='cancel')
        button_cancel.set_hexpand(True)
        button_cancel.connect("clicked", on_button_cancel)
        button_grid.add(button_cancel)

        def on_keypress(_, event):
            if event.keyval == Gdk.KEY_Escape:
                self.close()
            # this is bad when pressing return in sourceview
            #elif event.keyval == Gdk.KEY_Return:
            #    on_button_ok(_)
        self.connect('key-press-event', on_keypress)

        self.show_all()


class NodeReplaceWindow(Gtk.Window):
    def __init__(self, node, main_window):
        super(NodeReplaceWindow, self).__init__()

        grid = Gtk.Grid()
        self.add(grid)

        (possible_replacements_description,
         possible_replacements) = node.get_possible_replacements()

        description = Gtk.Label(label=possible_replacements_description)
        grid.add(description)

        if [r for r in possible_replacements if r.name == 'Integer']:
            # Integer
            grid_entry = Gtk.Grid()
            grid_entry.add(Gtk.Label(label='Integer: '))
            entry = Gtk.Entry()
            if not isinstance(node, PlaceholderNode):
                text = str(node.name)
            else:
                text = '0'
            entry.set_text(text)
            grid_entry.add(entry)
            grid.attach_next_to(grid_entry, description,
                                Gtk.PositionType.BOTTOM, 1, 1)
            # grid.add(grid_entry)

            def on_button_replace(_):
                integer = entry.get_text()
                try:
                    int(integer)
                except:
                    main_window.log_append_message(
                        Error(str(integer) + ' is not a valid Integer'))
                    return
                replacement = node.get_int_replacement(integer)
                ret = main_window.cpp_tree.replace_node(node, replacement)
                main_window.log_append_message(ret)
                main_window.redraw_all()
                self.close()

            button_replace = Gtk.Button(label='replace')
            button_replace.connect("clicked", on_button_replace)
            grid.attach_next_to(button_replace, grid_entry,
                                Gtk.PositionType.BOTTOM, 1, 1)

            def on_keypress(_, event):
                if event.keyval == Gdk.KEY_Escape:
                    self.close()
                elif event.keyval == Gdk.KEY_Return:
                    on_button_replace(_)
            self.connect('key-press-event', on_keypress)
        else:
            # no Integer
            def listbox_create_widget(node_line):
                listbox_row = ListBoxRowWithNode()
                listbox_row.add_node(node_line.node)
                grid = Gtk.Grid()
                grid.add(Gtk.Label(label=node_line.node.name))
                listbox_row.add(grid)
                return listbox_row

            store = Gio.ListStore()
            listbox = Gtk.ListBox()
            listbox.set_vexpand(True)
            listbox.set_hexpand(True)
            listbox.set_selection_mode(Gtk.SelectionMode.SINGLE)
            listbox.bind_model(store, listbox_create_widget)

            scroll = Gtk.ScrolledWindow()
            scroll.add(listbox)
            scroll.set_min_content_height(500)
            scroll.set_min_content_width(500)
            grid.attach_next_to(scroll, description,
                                Gtk.PositionType.BOTTOM, 1, 1)
            # grid.add(scroll)

            def on_button_replace(_):
                replacement = listbox.get_selected_row().node
                #replacement.add_missing_default_python_settings(
                #    main_window.cpp_tree.undo_stack.get_current_root().settings_dict)
                ret = main_window.cpp_tree.replace_node(node, replacement)
                main_window.log_append_message(ret)
                main_window.redraw_all()
                self.close()

            button_replace = Gtk.Button(label='replace')
            button_replace.connect("clicked", on_button_replace)

            listbox.set_activate_on_single_click(False)

            def row_double_clicked(_, row):
                on_button_replace(_)
            listbox.connect('row-activated', row_double_clicked)

            grid.attach_next_to(button_replace, scroll,
                                Gtk.PositionType.BOTTOM, 1, 1)

            for replacement in possible_replacements:
                store.append(NodeLine(replacement, 0))

            def on_keypress(_, event):
                if event.keyval == Gdk.KEY_Escape:
                    self.close()
                elif event.keyval == Gdk.KEY_Return:
                    on_button_replace(_)
            self.connect('key-press-event', on_keypress)

        self.show_all()


class MainWindow(Gtk.ApplicationWindow):
    def __init__(self, app):
        Gtk.Window.__init__(self, application=app)
        self.app = app
        self.init_ui()
        self.init_backend()

    def init_backend(self):
        self.cpp_tree = CPPTree()
        self.load_empty_simulation()

    def load_empty_simulation(self):
        ret = self.cpp_tree.load_empty_simulation()
        self.log_append_message(ret)
        self.redraw_all()

    def redraw_python(self):
        self.redraw_textview_python_code()
        self.redraw_treeview_python()

    def redraw_all(self):
        self.redraw_textview_cpp_code()
        self.redraw_treeview_cpp()
        self.redraw_python()

    def redraw_textview_python_code(self):
        text = str(self.cpp_tree.get_python_settings())
        self.text_view_python_code.get_buffer().set_text(text)

    def redraw_textview_cpp_code(self):
        text = str(self.cpp_tree)
        self.text_view_cpp_code.get_buffer().set_text(text)

    def redraw_treeview_cpp(self):
        # get currently selected index
        selected_row = self.cpp_treeview_listbox.get_selected_row()
        selected_index = 0
        while True:
            row = self.cpp_treeview_listbox.get_row_at_index(selected_index)
            if selected_row == row:
                break
            selected_index = selected_index + 1

        self.cpp_treeview_store.remove_all()
        self.redraw_treeview_cpp_recursive(
            self.cpp_tree.undo_stack.get_current_root(), 0)
        self.cpp_treeview_listbox.select_row(self.cpp_treeview_listbox.get_row_at_index(0))

        # reselect row at index
        row = self.cpp_treeview_listbox.get_row_at_index(selected_index)
        self.cpp_treeview_listbox.select_row(row)

    def redraw_treeview_cpp_recursive(self, node, depth):
        self.cpp_treeview_store.append(NodeLine(node, depth))
        for child in node.childs.get_childs():
            self.redraw_treeview_cpp_recursive(child, depth + 1)

    def redraw_treeview_python(self):
        self.redraw_treeview_python_code()
        self.redraw_treeview_python_list()

    def redraw_treeview_python_code(self):
        try:
            row = self.cpp_treeview_listbox.get_selected_row()
            self.cpp_treeview_current_row = row
            node = row.node
            text = str(node.settings_dict)
        except:
            text = ''
            self.cpp_treeview_current_row = None
        self.python_treeview_code.get_buffer().set_text(text)

    def redraw_treeview_python_list(self):
        self.python_treeview_store.remove_all()
        try:
            node = self.cpp_treeview_listbox.get_selected_row().node
            self.redraw_treeview_python_list_recursive(node.settings_dict, 0)
            self.python_treeview_listbox.select_row(self.python_treeview_listbox.get_row_at_index(0))
        except: pass

    def redraw_treeview_python_list_recursive(self, settings, depth, parent=None):
        if isinstance(settings, SettingsList) or isinstance(settings, SettingsDict):
            for e in settings:
                self.redraw_treeview_python_list_recursive(e, depth, parent=settings)
        elif isinstance(settings, SettingsDictEntry) or isinstance(settings, SettingsListEntry):
            self.python_treeview_store.append(SettingsLine(settings, depth, parent=parent))
            if not isinstance(settings.value, str):
                self.redraw_treeview_python_list_recursive(settings.value, depth + 1, parent=settings)
        elif isinstance(settings, SettingsChildPlaceholder):
            self.python_treeview_store.append(SettingsLine(settings, depth, parent=parent))

    def open_file(self, type):
        if type == 'cpp':
            title = "Please choose a c++ file"
        else:
            title = "Please choose a python file"
        dialog = Gtk.FileChooserDialog(
            title=title, parent=self, action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(
            Gtk.STOCK_CANCEL,
            Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN,
            Gtk.ResponseType.OK,
        )

        filter_type = Gtk.FileFilter()
        if type == 'cpp':
            filter_type.set_name("C++ files")
            filter_type.add_mime_type("text/x-c")
        else:
            filter_type.set_name("Python files")
            filter_type.add_mime_type("text/x-python")
        dialog.add_filter(filter_type)

        filter_any = Gtk.FileFilter()
        filter_any.set_name("Any files")
        filter_any.add_pattern("*")
        dialog.add_filter(filter_any)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            try:
                file = open(dialog.get_filename(), "r")
                text = file.read()
                file.close()
                if type == 'cpp':
                    self.text_view_cpp_code.get_buffer().set_text(text)
                    self.on_button_apply_cpp_code(None)
                else:
                    self.text_view_python_code.get_buffer().set_text(text)
                    self.on_button_apply_python_code(None)
                self.log_append_message(
                    Info('loaded ' + dialog.get_filename()))
            except:
                pass
        dialog.destroy()

    def save_file(self, type):
        if type == 'cpp':
            title = "Save to c++ file"
        else:
            title = "Save to python file"
        dialog = Gtk.FileChooserDialog(
            title=title, parent=self, action=Gtk.FileChooserAction.SAVE)
        dialog.set_do_overwrite_confirmation(True)
        dialog.add_buttons(
            Gtk.STOCK_CANCEL,
            Gtk.ResponseType.CANCEL,
            Gtk.STOCK_SAVE,
            Gtk.ResponseType.OK,
        )

        filter_type = Gtk.FileFilter()
        if type == 'cpp':
            filter_type.set_name("C++ files")
            filter_type.add_mime_type("text/x-c")
        else:
            filter_type.set_name("Python files")
            filter_type.add_mime_type("text/x-python")
        dialog.add_filter(filter_type)

        filter_any = Gtk.FileFilter()
        filter_any.set_name("Any files")
        filter_any.add_pattern("*")
        dialog.add_filter(filter_any)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            try:
                file = open(dialog.get_filename(), "w")
                if type == 'cpp':
                    text_bounds = self.text_view_cpp_code.get_buffer().get_bounds()
                    text = self.text_view_cpp_code.get_buffer().get_text(
                        text_bounds[0], text_bounds[1], True)
                    file.write(text)
                    self.log_append_message(
                        Info('saved c++-code to ' + dialog.get_filename()))
                else:
                    text_bounds = self.text_view_python_code.get_buffer().get_bounds()
                    text = self.text_view_python_code.get_buffer().get_text(
                        text_bounds[0], text_bounds[1], True)
                    file.write(text)
                    self.log_append_message(
                        Info('saved python-code to ' + dialog.get_filename()))
                file.close()
            except:
                pass
        dialog.destroy()

    # binded to self.cpp_treeview_listbox
    # if we append item to self.cpp_treeview_store this function gets called and creates the widget

    def cpp_treeview_listbox_create_widget(self, list_node):
        node = list_node.node
        depth = list_node.depth

        grid = Gtk.Grid()
        for _ in range(depth):
            grid.add(Gtk.Label(label='  '))
            #grid.add(Gtk.Label(label='   '))

        description = node.get_contextual_description()
        if not isinstance(node, RootNode):
            label = Gtk.Label(label=description + ' ')
            grid.add(label)

        if isinstance(node, PlaceholderNode):
            if node.needed:
                label = 'add necessary template'
                color = Gdk.RGBA(1, 0, 0, .5)
                grid.override_background_color(Gtk.StateType.NORMAL, color)
            else:
                label = 'add optional template'

            def on_button_add_node(_):
                self.cpp_treeview_replace_node(node)

            button_add_node = Gtk.Button(label=label)
            button_add_node.connect("clicked", on_button_add_node)
            grid.add(button_add_node)
        else:
            label = Gtk.Label()
            label.set_markup("<b>" + str(node.name) + "</b>")
            grid.add(label)
            if not isinstance(node, RootNode):
                spacer = Gtk.Label(label=' ')
                grid.add(spacer)
                button_delete_node = Gtk.Button()
                icon = Gio.ThemedIcon(name="edit-delete-symbolic")
                image = Gtk.Image.new_from_gicon(icon, Gtk.IconSize.BUTTON)
                button_delete_node.add(image)

                def on_button_delete_node(_):
                    self.cpp_treeview_delete_node(node)
                button_delete_node.connect("clicked", on_button_delete_node)
                grid.add(button_delete_node)

            errors = node.validate_cpp_src()
            if errors:
                s = ''
                for error in errors:
                    s = s + str(error) + '\n'
                grid.set_tooltip_text(s[:-1])

                color = Gdk.RGBA(1, 0, 0, .5)
                grid.override_background_color(Gtk.StateType.NORMAL, color)
        grid.show_all()
        row = ListBoxRowWithNode()
        row.node = node
        row.add(grid)
        return row

    def cpp_treeview_replace_node(self, node):
        if self.check_python_treeview_for_unapplied_code():
            return
        _window = NodeReplaceWindow(node, self)

    def cpp_treeview_delete_node(self, node):
        if self.cpp_treeview_current_row != None and node != self.cpp_treeview_current_row.node:
            if self.check_python_treeview_for_unapplied_code():
                return
        ret = self.cpp_tree.delete_node(node)
        self.log_append_message(ret)
        self.redraw_all()

    def python_treeview_listbox_create_widget(self, list_settings):
        settings = list_settings.settings
        depth = list_settings.depth
        parent = list_settings.parent
        grid = Gtk.Grid()
        for _ in range(depth):
            #grid.add(Gtk.Label(label='   '))
            grid.add(Gtk.Label(label='  '))
        try:
            if settings.comments:
                grid.set_tooltip_text('comments: ' + str(settings.comments) + '\ndefault-comment: ' + str(settings.default_comment))
            else:
                grid.set_tooltip_text('default-comment: ' + str(settings.default_comment))
        except: pass
        activated = True
        try:
            if settings.activated == False:
                activated = False
        except: pass
        if isinstance(settings, SettingsDict):
            grid.add(Gtk.Label(label='dict:'))
        if isinstance(settings, SettingsList):
            grid.add(Gtk.Label(label='list:'))
        elif isinstance(settings, SettingsDictEntry):
            if activated == False:
                # set textcolor to grey if this entry is not activated
                grid.override_color(Gtk.StateFlags.NORMAL, Gdk.RGBA(0.5, 0.5, 0.5, 1.0))

            checkbox_activated = Gtk.CheckButton()
            # with the callback function toggled_cb
            if activated:
                checkbox_activated.set_active(True)
            def toggled_checkbox(button):
                self.cpp_tree.undo_stack.duplicate_current_state()
                if button.get_active():
                    # activate this and parents
                    settings.activate_recursive()
                else:
                    settings.activated = False
                    pass
                self.redraw_python()
            checkbox_activated.connect("toggled", toggled_checkbox)
            grid.add(checkbox_activated)

            v = ''
            if isinstance(settings.value, str):
                v = settings.value
            elif isinstance(settings.value, SettingsDict):
                v = '{..}'
            elif isinstance(settings.value, SettingsList):
                v = '[..]'
            grid.add(Gtk.Label(label=settings.key + ' : ' + v))
            if settings.is_unknown:
                color = Gdk.RGBA(1, 1, 0, .5)
                grid.override_background_color(Gtk.StateType.NORMAL, color)
                grid.set_tooltip_text(settings.key + ' is an unknown setting in this context')

        elif isinstance(settings, SettingsListEntry):
            v = ''
            if isinstance(settings.value, str):
                v = settings.value
            elif isinstance(settings.value, SettingsDict):
                v = '{..}'
            elif isinstance(settings.value, SettingsList):
                v = '[..]'
            grid.add(Gtk.Label(label=v))
        elif isinstance(settings, SettingsChildPlaceholder):
            grid.add(Gtk.Label(label=settings.comment))

        try:
            doc_link = settings.doc_link
            if not webmode and not doc_link == None:
                grid.add(Gtk.Label(label=' '))
                button_open_doc = Gtk.Button()
                button_open_doc.set_tooltip_text('open web-documentation')
                #icon = Gio.ThemedIcon(name="help-about")
                icon = Gio.ThemedIcon(name="help-faq")
                #icon = Gio.ThemedIcon(name="dialog-question")
                image = Gtk.Image.new_from_gicon(icon, Gtk.IconSize.BUTTON)
                button_open_doc.add(image)
                def on_button_open_doc(_):
                    if sys.platform=='darwin':
                        subprocess.Popen(['open', doc_link])
                    elif sys.platform == 'win32':
                        os.startfile(doc_link)
                    else:
                        subprocess.Popen(['xdg-open', doc_link])
                button_open_doc.connect("clicked", on_button_open_doc)
                grid.add(button_open_doc)
        except: pass

        if activated and not isinstance(settings, SettingsChildPlaceholder):
            grid.add(Gtk.Label(label=' '))
            button_delete_settings = Gtk.Button()
            icon = Gio.ThemedIcon(name="edit-delete-symbolic")
            image = Gtk.Image.new_from_gicon(icon, Gtk.IconSize.BUTTON)
            button_delete_settings.add(image)
            def on_button_delete_settings(_):
                if self.check_python_treeview_for_unapplied_code():
                    return
                self.cpp_tree.undo_stack.duplicate_current_state()
                parent.remove(settings)
                self.redraw_python()
                node = self.cpp_treeview_listbox.get_selected_row().node
                text_bounds = self.python_treeview_code.get_buffer().get_bounds()
                text = self.python_treeview_code.get_buffer().get_text(
                    text_bounds[0], text_bounds[1], True)
                self.cpp_tree.parse_python_settings(text, node, undoable=False)
                self.redraw_python()
                self.log_append_message(Info('removed ' + GLib.markup_escape_text(str(type(settings))) + ' from python-options'))
            button_delete_settings.connect("clicked", on_button_delete_settings)
            grid.add(button_delete_settings)

        grid.show_all()
        row = ListBoxRowWithSettings()
        row.settings = settings
        row.add(grid)
        return row

    # returns True if there are unapplied changes in python_treeview_code, that should not be discarded
    # returns Fals otherwise
    def check_python_treeview_for_unapplied_code(self):
        text_bounds = self.python_treeview_code.get_buffer().get_bounds()
        text = self.python_treeview_code.get_buffer().get_text(
            text_bounds[0], text_bounds[1], True)
        if self.cpp_treeview_current_row != None and str(self.cpp_treeview_current_row.node.settings_dict) != text:
            dialog = DiscardNodeChangesDialog(self)
            response = dialog.run()
            dialog.destroy()
            if response == Gtk.ResponseType.CANCEL or response == Gtk.ResponseType.DELETE_EVENT:
                self.cpp_treeview_listbox.select_row(
                    self.cpp_treeview_current_row)
                return True
        return False

    def log_append_message(self, message):
        if isinstance(message, list):
            for m in message:
                self.log_append_line(str(m), m.color)
        else:
            self.log_append_line(str(message), message.color)

    def log_append_line(self, text, color=None):
        buffer = self.text_view_log.get_buffer()
        iter = buffer.get_end_iter()
        col = ''
        if color:
            col = ' color="' + color + '"'
        for line in str(text).splitlines():
            if buffer.get_char_count() > 0:
                buffer.insert(iter, '\n')
            # move mark
            buffer.move_mark(self.log_text_mark_end, iter)
            buffer.insert_markup(iter, '<span' + col +
                                 '>' + line + '</span>', -1)
        # scroll to end
        self.text_view_log.scroll_to_mark(
            self.log_text_mark_end, 0, False, 0, 0)

    def on_button_add_defaults_python_code(self, _):
        if self.check_python_treeview_for_unapplied_code():
            return
        ret = self.cpp_tree.activate_all_default_python_settings()
        self.log_append_message(ret)
        self.redraw_all()

    def on_button_apply_python_code(self, _):
        if self.check_python_treeview_for_unapplied_code():
            return
        text_bounds = self.text_view_python_code.get_buffer().get_bounds()
        text = self.text_view_python_code.get_buffer().get_text(
            text_bounds[0], text_bounds[1], True)
        rets = self.cpp_tree.parse_python_settings(text)
        self.log_append_message(rets)
        #if not isinstance(rets, Error):
        #    self.redraw_python()
        self.redraw_all()

    def on_button_apply_cpp_code(self, _):
        if self.check_python_treeview_for_unapplied_code():
            return
        text_bounds = self.text_view_cpp_code.get_buffer().get_bounds()
        text = self.text_view_cpp_code.get_buffer().get_text(
            text_bounds[0], text_bounds[1], True)
        rets = self.cpp_tree.parse_cpp_src(
            text, validate_semantics=self.checkbox_validate_semantics.get_active())
        self.log_append_message(rets)
        if not any(isinstance(ret, Error) or isinstance(ret, Warning) for ret in rets):
            self.log_append_message(
                Info('trying to parse python-settings from the old cpp-tree to the new one'))
            text_bounds = self.text_view_python_code.get_buffer().get_bounds()
            text = self.text_view_python_code.get_buffer().get_text(
                text_bounds[0], text_bounds[1], True)
            rets = self.cpp_tree.parse_python_settings(text)
            self.log_append_message(rets)
            self.redraw_all()

    def on_button_undo(self, _):
        ret = self.cpp_tree.undo_stack.undo()
        self.log_append_message(ret)
        if not isinstance(ret, Error):
            self.redraw_all()

    def on_button_redo(self, _):
        ret = self.cpp_tree.undo_stack.redo()
        self.log_append_message(ret)
        if not isinstance(ret, Error):
            self.redraw_all()

    # def on_button_reset(self, _):
    #    self.load_empty_simulation()

    def on_python_treeview_button_apply(self, _):
        try:
            node = self.cpp_treeview_listbox.get_selected_row().node

            text_bounds = self.python_treeview_code.get_buffer().get_bounds()
            text = self.python_treeview_code.get_buffer().get_text(
                text_bounds[0], text_bounds[1], True)

            rets = self.cpp_tree.parse_python_settings(text, node)
            self.log_append_message(rets)
            if not isinstance(rets, Error):
                self.redraw_python()
            else:
                self.redraw_all()
        except:
            self.log_append_message(
                Error('Can\'t apply settings if no Node is selected'))

    def on_python_treeview_button_add_defaults(self, _):
        if self.check_python_treeview_for_unapplied_code():
            return
        try:
            node = self.cpp_treeview_listbox.get_selected_row().node

            rets = self.cpp_tree.activate_all_default_python_settings(node)
            self.log_append_message(rets)
            self.redraw_all()
        except:
            self.log_append_message(
                Error('Can\'t add default settings if no Node is selected'))

    #def on_python_treeview_button_add_missing_default_setting(self, _):
    #    if self.check_python_treeview_for_unapplied_code():
    #        return
    #    #try:
    #    if 1 == 1:
    #        node = self.cpp_treeview_listbox.get_selected_row().node
    #        missing_settings = node.get_unused_python_settings()
    #        print(missing_settings)
    #        if missing_settings:
    #            pass
    #            # TODO create window with a list of missing_settings and let the user select one
    #            # TODO then add the selected setting (and needed parents to self.settings_dict)
    #            # TODO don't forget to create undo
    #            # TODO also check if a key is in SettingsConditionals while adding with has_key() and get_value()
    #        else:
    #            self.log_append_message(Info('There is no additional setting, that you could add'))
    #    #except:
    #    #    self.log_append_message(
    #    #        Error('Can\'t add missing setting if no Node is selected'))

    def init_ui(self):
        self.cpp_treeview_current_row = None

        self.set_title("opendihu - webapp")
        self.connect("destroy", Gtk.main_quit)

        # menu_bar
        action = Gio.SimpleAction.new("menu_new", None)

        def on_new(_action, _parameter):
            self.load_empty_simulation()
        action.connect("activate", on_new)
        self.add_action(action)

        action = Gio.SimpleAction.new("menu_undo", None)

        def on_undo(_action, _parameter):
            self.on_button_undo(None)
        action.connect("activate", on_undo)
        self.add_action(action)

        action = Gio.SimpleAction.new("menu_redo", None)

        def on_redo(_action, _parameter):
            self.on_button_redo(None)
        action.connect("activate", on_redo)
        self.add_action(action)

        if not webmode:
            action = Gio.SimpleAction.new("menu_open_cpp", None)

            def on_open_cpp(_action, _parameter):
                self.open_file('cpp')
            action.connect("activate", on_open_cpp)
            self.add_action(action)

            action = Gio.SimpleAction.new("menu_open_python", None)

            def on_open_python(_action, _parameter):
                self.open_file('python')
            action.connect("activate", on_open_python)
            self.add_action(action)

            action = Gio.SimpleAction.new("menu_save_cpp", None)

            def on_save_cpp(_action, _parameter):
                self.save_file('cpp')
            action.connect("activate", on_save_cpp)
            self.add_action(action)

            action = Gio.SimpleAction.new("menu_save_python", None)

            def on_save_python(_action, _parameter):
                self.save_file('python')
            action.connect("activate", on_save_python)
            self.add_action(action)

        # header_bar
        self.header_bar = Gtk.HeaderBar()
        self.set_titlebar(self.header_bar)

        self.button_undo = Gtk.Button()
        icon = Gio.ThemedIcon(name="edit-undo-symbolic")
        image = Gtk.Image.new_from_gicon(icon, Gtk.IconSize.BUTTON)
        self.button_undo.add(image)
        self.button_undo.connect("clicked", self.on_button_undo)
        self.header_bar.pack_start(self.button_undo)

        self.button_redo = Gtk.Button()
        icon = Gio.ThemedIcon(name="edit-redo-symbolic")
        image = Gtk.Image.new_from_gicon(icon, Gtk.IconSize.BUTTON)
        self.button_redo.add(image)
        self.button_redo.connect("clicked", self.on_button_redo)
        self.header_bar.pack_start(self.button_redo)

        #self.button_reset = Gtk.Button(label='load empty simulation')
        #self.button_reset.connect("clicked", self.on_button_reset)
        # self.header_bar.pack_start(self.button_reset)

        # main grid
        self.paned_main = Gtk.Paned.new(Gtk.Orientation.VERTICAL)
        self.add(self.paned_main)

        # log (lower)
        self.text_view_log = Gtk.TextView()
        self.text_view_log.set_editable(False)
        buffer = self.text_view_log.get_buffer()
        iter = buffer.get_end_iter()
        self.log_text_mark_end = buffer.create_mark('the-end', iter, True)
        self.scroll_log = Gtk.ScrolledWindow()
        self.scroll_log.set_min_content_height(100)
        self.scroll_log.add(self.text_view_log)
        self.tabs_log = Gtk.Notebook()
        self.tabs_log.append_page(self.scroll_log, Gtk.Label(label='log'))
        self.paned_main.pack2(self.tabs_log, resize=False, shrink=False)

        # upper tabs
        self.tabs_main = Gtk.Notebook()
        self.paned_main.pack1(self.tabs_main, resize=True, shrink=False)

        # code view
        self.grid_code_view = Gtk.Grid(column_homogeneous=True)
        self.tabs_main.append_page(
            self.grid_code_view, Gtk.Label(label='Global-View'))

        # cpp code view
        self.tabs_cpp_code = Gtk.Notebook()
        self.grid_code_view.add(self.tabs_cpp_code)

        self.grid_cpp_code = Gtk.Grid()

        self.text_view_cpp_code = GtkSource.View()
        language_manager = GtkSource.LanguageManager()
        self.text_view_cpp_code.get_buffer().set_language(
            language_manager.get_language('cpp'))
        self.text_view_cpp_code.set_vexpand(True)
        self.text_view_cpp_code.set_hexpand(True)
        self.scroll_cpp_code = Gtk.ScrolledWindow()
        self.scroll_cpp_code.add(self.text_view_cpp_code)
        self.grid_cpp_code.add(self.scroll_cpp_code)

        self.grid_cpp_code_buttons = Gtk.Grid()
        self.button_apply_cpp_code = Gtk.Button(label='apply changes')
        self.button_apply_cpp_code.connect(
            "clicked", self.on_button_apply_cpp_code)
        #self.button_verify_cpp_code = Gtk.Button(label='validate')
        #self.button_verify_cpp_code.connect("clicked", self.on_button_verify_cpp_code)
        self.checkbox_validate_semantics = Gtk.CheckButton()
        self.checkbox_validate_semantics.set_label(
            "check code for semantic errors")
        self.checkbox_validate_semantics.set_active(True)
        self.grid_cpp_code_buttons.add(self.button_apply_cpp_code)
        self.grid_cpp_code_buttons.attach_next_to(
            self.checkbox_validate_semantics, self.button_apply_cpp_code, Gtk.PositionType.RIGHT, 1, 1)
        self.grid_cpp_code.attach_next_to(
            self.grid_cpp_code_buttons, self.scroll_cpp_code, Gtk.PositionType.BOTTOM, 1, 1)

        self.tabs_cpp_code.append_page(
            self.grid_cpp_code, Gtk.Label(label='C++ Code'))

        # python code view
        self.tabs_python_code = Gtk.Notebook()
        self.grid_code_view.attach_next_to(
            self.tabs_python_code, self.tabs_cpp_code, Gtk.PositionType.RIGHT, 1, 1)

        self.grid_python_code = Gtk.Grid()
        self.tabs_python_code.append_page(
            self.grid_python_code, Gtk.Label(label='Global Python-Code'))

        self.text_view_python_code = GtkSource.View()
        self.text_view_python_code.get_buffer().set_language(
            language_manager.get_language('python3'))
        self.text_view_python_code.set_vexpand(True)
        self.text_view_python_code.set_hexpand(True)
        self.scroll_python_code = Gtk.ScrolledWindow()
        self.scroll_python_code.add(self.text_view_python_code)
        self.grid_python_code.add(self.scroll_python_code)

        self.grid_python_code_buttons = Gtk.Grid()
        self.grid_python_code.attach_next_to(
            self.grid_python_code_buttons, self.scroll_python_code, Gtk.PositionType.BOTTOM, 1, 1)

        self.button_apply_python_code = Gtk.Button(label='apply changes')
        self.button_apply_python_code.connect(
            "clicked", self.on_button_apply_python_code)
        self.grid_python_code_buttons.add(self.button_apply_python_code)

        self.button_add_defaults_python_code = Gtk.Button(
            label='activate all possible default settings recursively')
        self.button_add_defaults_python_code.connect(
            "clicked", self.on_button_add_defaults_python_code)
        self.grid_python_code_buttons.attach_next_to(
            self.button_add_defaults_python_code, self.button_apply_python_code, Gtk.PositionType.RIGHT, 1, 1)

        # tree view
        self.grid_treeview = Gtk.Grid(column_homogeneous=True)
        self.tabs_main.prepend_page(self.grid_treeview, Gtk.Label(label='Tree-View'))

        # cpp tree view
        self.tabs_cpp_treeview = Gtk.Notebook()
        self.grid_treeview.add(self.tabs_cpp_treeview)

        self.cpp_treeview_store = Gio.ListStore()
        self.cpp_treeview_listbox = Gtk.ListBox()
        self.cpp_treeview_listbox.set_vexpand(True)
        self.cpp_treeview_listbox.set_hexpand(True)
        self.cpp_treeview_listbox.set_selection_mode(Gtk.SelectionMode.SINGLE)
        self.cpp_treeview_listbox.bind_model(
            self.cpp_treeview_store, self.cpp_treeview_listbox_create_widget)

        self.cpp_treeview_listbox.set_activate_on_single_click(False)

        def row_double_clicked(_, row):
            node = self.cpp_treeview_listbox.get_selected_row().node
            self.cpp_treeview_replace_node(node)
        self.cpp_treeview_listbox.connect('row-activated', row_double_clicked)

        def row_clicked(_, row):
            if row != None and self.cpp_treeview_current_row != None:
                if row.node != self.cpp_treeview_current_row.node:
                    # check if we have changes that are not applied yet
                    if self.check_python_treeview_for_unapplied_code():
                        self.cpp_treeview_listbox.select_row(
                            self.cpp_treeview_current_row)
                        return
                else:
                    return
            self.redraw_treeview_python()
        self.cpp_treeview_listbox.connect('row-selected', row_clicked)

        self.scroll_cpp_treeview = Gtk.ScrolledWindow()
        self.scroll_cpp_treeview.add(self.cpp_treeview_listbox)
        self.tabs_cpp_treeview.append_page(
            self.scroll_cpp_treeview, Gtk.Label(label='C++ Tree'))


        # python treeview list-tab
        self.tabs_python_treeview = Gtk.Notebook()
        self.python_treeview_store = Gio.ListStore()
        self.python_treeview_listbox = Gtk.ListBox()
        self.python_treeview_listbox.set_vexpand(True)
        self.python_treeview_listbox.set_hexpand(True)
        self.python_treeview_listbox.set_selection_mode(Gtk.SelectionMode.SINGLE)
        self.python_treeview_listbox.bind_model(
            self.python_treeview_store, self.python_treeview_listbox_create_widget)

        self.python_treeview_listbox.set_activate_on_single_click(False)

        def python_row_double_clicked(_, row):
            if self.check_python_treeview_for_unapplied_code():
                return
            _window = PythonSettingsChangeWindow(row.settings, self)

        self.python_treeview_listbox.connect('row-activated', python_row_double_clicked)

        def python_row_clicked(_, row):
            pass
        self.python_treeview_listbox.connect('row-selected', python_row_clicked)

        self.scroll_python_treeview = Gtk.ScrolledWindow()
        self.scroll_python_treeview.add(self.python_treeview_listbox)
        self.tabs_python_treeview.append_page(
            self.scroll_python_treeview, Gtk.Label(label='Selected Python-Tree'))


        # python treeview code
        self.python_treeview_outer_grid = Gtk.Grid()
        self.python_treeview_outer_grid.add(self.tabs_python_treeview)
        self.grid_treeview.attach_next_to(
            self.python_treeview_outer_grid, self.tabs_cpp_treeview, Gtk.PositionType.RIGHT, 1, 1)
        self.python_treeview_grid = Gtk.Grid()
        self.tabs_python_treeview.append_page(
            self.python_treeview_grid, Gtk.Label(label='Selected Python-Code'))
        self.python_treeview_scroll = Gtk.ScrolledWindow()
        self.python_treeview_grid.add(self.python_treeview_scroll)

        self.python_treeview_code = GtkSource.View()
        self.python_treeview_code.get_buffer().set_language(
            language_manager.get_language('python3'))
        self.python_treeview_code.set_vexpand(True)
        self.python_treeview_code.set_hexpand(True)
        self.python_treeview_scroll.add(self.python_treeview_code)

        self.python_treeview_buttons = Gtk.Grid()
        self.python_treeview_outer_grid.attach_next_to(
            self.python_treeview_buttons, self.tabs_python_treeview, Gtk.PositionType.BOTTOM, 1, 1)

        self.python_treeview_code_buttons = Gtk.Grid()
        self.python_treeview_grid.attach_next_to(
            self.python_treeview_code_buttons, self.python_treeview_scroll, Gtk.PositionType.BOTTOM, 1, 1)

        self.python_treeview_button_apply = Gtk.Button(label='apply changes')
        self.python_treeview_button_apply.connect(
            "clicked", self.on_python_treeview_button_apply)
        self.python_treeview_code_buttons.add(self.python_treeview_button_apply)

        self.python_treeview_button_add_defaults = Gtk.Button(
            label='activate all possible default settings on selected node')
        self.python_treeview_button_add_defaults.connect(
            "clicked", self.on_python_treeview_button_add_defaults)
        self.python_treeview_buttons.add(
            self.python_treeview_button_add_defaults)



class MainApplication(Gtk.Application):
    def __init__(self):
        Gtk.Application.__init__(self)

    # override super.run to process sys.argv
    # TODO can gtk handle argv for us instead?
    def run(self, argv):
        try:
            if str(sys.argv[1]) == '--web':
                # set webmode if first argument is --web
                global webmode
                webmode = True
        except: pass
        return super().run([])

    def do_activate(self):
        win = MainWindow(self)
        win.show_all()
        win.maximize()

    def do_startup(self):
        Gtk.Application.do_startup(self)

        # menu_bar
        builder = Gtk.Builder()
        builder.add_from_file("menubar.ui")
        self.set_menubar(builder.get_object("menubar"))

        if not webmode:
            action = Gio.SimpleAction.new("menu_quit", None)

            def on_quit(_action, _parameter):
                sys.exit()
            action.connect("activate", on_quit)
            self.add_action(action)


webmode = False
app = MainApplication()
exit_status = app.run(sys.argv)
sys.exit(exit_status)
