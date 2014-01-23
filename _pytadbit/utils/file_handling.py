"""
22 Jan 2014


"""

def check_pik(path):
    with open(path, "r") as f:
        f.seek (0, 2)           # Seek @ EOF
        fsize = f.tell()        # Get Size
        f.seek (max (fsize-2, 0), 0) # Set pos @ last n chars
        key = f.read()       # Read to end
    return key == 's.'
