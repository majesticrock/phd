import string
def alpha_suffix(array):
    # Příprava seznamu sufixů: začínáme s prázdným řetězcem, pak přidáme "_a", "_b", ...
    suffixes = [""] + [f"_{char}" for char in string.ascii_lowercase]
    for i, item in enumerate(array):
        yield (item, suffixes[i])
        
def alpha_label(axes, array):
    label = [f"({char})" for char in string.ascii_lowercase]
    for i, item in enumerate(array):
        yield (axes[i], item, label[i])