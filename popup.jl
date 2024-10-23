using Gtk

function new_resdir_gui()
    win = GtkWindow("New result directory",300,50)

    ent = GtkEntry()
    set_gtk_property!(ent,:text,"Type full directory: ")

    push!(win,ent)
    Gtk.showall(win)


    str = ""
    if isinteractive()
        c = Condition()
        signal_connect(win, "key-press-event") do widget, event
            grab_focus(widget)
            if event.keyval == 65293 #Enter
                str = get_gtk_property(ent,:text,String)
                notify(c)
            end
        end
        wait(c)
    end
    present(win) 

    Gtk.destroy(win)

    return str
end

function dir_not_empty()
    win = GtkWindow("resdir is not empty. Define new dir?", 400, 200)

    a = GtkButton("No")
    push!(win,a)

    #b = GtkButton("Yes")
    #push!(win,b)

    showall(win)
end
