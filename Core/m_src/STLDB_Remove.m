function STLDB_Remove(id)
    global BreachGlobOpt
    if isKey(BreachGlobOpt.STLDB, id)
        BreachGlobOpt.STLDB.remove(id);
    end;