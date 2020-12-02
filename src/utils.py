def is_variable_constraint(lhs):
    return list(lhs.values()).count(0) == (len(lhs) - 1)

def get_var_from_constraint(lhs):
    return next(k for k, v in lhs.items() if v != 0)

def reorder_by_index(dict_, order):
    return {k: v for k, v in sorted(dict_.items(), key=lambda tup: order.index(tup[0]))}