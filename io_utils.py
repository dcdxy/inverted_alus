#!/usr/bin/env python

import os
import pickle
from base64 import b64decode
from io import BytesIO
import matplotlib.pyplot as plt
import PIL
from IPython import get_ipython
from IPython.core import magic_arguments
from IPython.core.magic import (Magics, cell_magic, magics_class)
from IPython.display import display
from IPython.utils.capture import capture_output

def pickle_enrich_dct(enrich_dct, pickle_filename, pickle_dir):
    """
    Ref: https://stackoverflow.com/questions/11218477
    """
    with open(os.path.join(pickle_dir, pickle_filename), 'wb') as handle:
        pickle.dump(enrich_dct, handle, protocol=pickle.HIGHEST_PROTOCOL)


def open_pickle_enrich_dct(pickle_filename, pickle_dir):
    """
    Ref: https://stackoverflow.com/questions/11218477
    """
    with open(os.path.join(pickle_dir, pickle_filename), 'rb') as handle:
        return pickle.load(handle)


# Ref: https://discourse.jupyter.org/t/cell-magic-to-save-image-output-as-a-png-file/11906/5
# Author: kolibril13
@magics_class
class MyMagic(Magics):

    @cell_magic
    @magic_arguments.magic_arguments()
    @magic_arguments.argument(
        "--path",
        "-p",
        default=None,
        help=("The path where the image will be saved to. When there is more then one image, multiple paths have to be defined"),
    )
    @magic_arguments.argument(
        "--compression",
        "-c",
        default=None,
        help=("Defines the amout of compression,  quality can be from 0.1 - 100 , images must be .jpg"),
    )
    def polaroid_camera(self, line, cell):
        args = magic_arguments.parse_argstring(MyMagic.polaroid_camera, line)
        paths = args.path.strip('"').split(' ')
        with capture_output(stdout=False, stderr=False, display=True) as result:
            self.shell.run_cell(cell) # thanks @krassowski for this idea!

        for output in result.outputs:
            display(output)
            data = output.data
            if 'image/png' in data:
                path = paths.pop(0)
                if not path:
                    raise ValueError('Too few paths given!')
                png_bytes = data['image/png']
                if isinstance(png_bytes, str):
                    png_bytes = b64decode(png_bytes)
                assert isinstance(png_bytes, bytes)
                bytes_io = BytesIO(png_bytes)
                img = PIL.Image.open(bytes_io)
                if args.compression:
                    if img.mode != "RGB":
                        img = img.convert("RGB")
                    img.save(path, "JPEG", optimize=True,
                             quality=int(args.compression))
                else:
                    img.save(path, 'png')


ipy = get_ipython()
ipy.register_magics(MyMagic)