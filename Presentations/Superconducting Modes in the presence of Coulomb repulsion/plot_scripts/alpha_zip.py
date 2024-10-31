import string
def alpha_suffix(array):
    # Příprava seznamu sufixů: začínáme s prázdným řetězcem, pak přidáme "_a", "_b", ...
    suffixes = [""] + [f"_{char}" for char in string.ascii_lowercase]
    for i, item in enumerate(array):
        yield (item, suffixes[i])
        
def alpha_label(axes, array, suffix="", begin=0):
    label = [f"({char}{suffix})" for char in string.ascii_lowercase]
    for i, item in enumerate(array):
        yield (axes[begin + i], item, label[begin + i])
