
def fetch_properties(method_class, ignore=[]):
    return [
        func
        for func in method_class.__dir__()
        if (not func.startswith("_"))
        & (hasattr(method_class, func))
        & (not func in ignore)
    ]