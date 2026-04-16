import os
import urllib.request

from RBGeneratorCpp import *
from rb_generator._version import __version__


def load_table(path=".", from_="github"):
    """
    Load rb22 and rb44 table data.

    If rb22.dat and rb44.dat already exist in `path`, they are used directly.
    Otherwise, data is acquired according to `from_`:
      - "github": download from GitHub releases
      - "local":  generate tables in-place using C++ generators

    Parameters
    ----------
    path : str
        Directory containing (or where to place) rb22.dat and rb44.dat.
    from_ : str
        Data source when files are not already present. "github" or "local".
    """
    rb22_path = os.path.join(path, "rb22.dat")
    rb44_path = os.path.join(path, "rb44.dat")

    if os.path.exists(rb22_path) and os.path.exists(rb44_path):
        # Files present — load directly
        rb22 = RB22()
        rb22.load_from_file(rb22_path)
        rb44 = RB44()
        rb44.load_from_file(rb44_path)
        return rb22, rb44

    os.makedirs(path, exist_ok=True)

    if from_ == "github":
        _download_from_github(path)
    elif from_ == "local":
        _generate_locally(path)
    else:
        raise ValueError(f"load_table: from_ must be 'github' or 'local', got {from_!r}")

    # Load after acquisition
    rb22 = RB22()
    rb22.load_from_file(rb22_path)
    rb44 = RB44()
    rb44.load_from_file(rb44_path)
    return rb22, rb44


def _download_from_github(path):
    url_base = (
        f"https://github.com/Agony5757/rb_generator"
        f"/releases/download/{__version__}"
    )
    for name in ("rb22.dat", "rb44.dat"):
        dest = os.path.join(path, name)
        if not os.path.exists(dest):
            url = f"{url_base}/{name}"
            print(f"Downloading {url} ...")
            urllib.request.urlretrieve(url, dest)


def _generate_locally(path):
    print("Generating rb22.dat ...")
    generate_table22(os.path.join(path, "rb22.dat"))
    print("Generating rb44.dat ...")
    generate_table44(os.path.join(path, "rb44.dat"))


__all__ = ["__version__", "load_table"]
