#!/usr/bin/env python3

import subprocess
import argparse
import PIL.Image
from PIL import Image, ImageChops, ImageDraw


# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("T", type=int)
ap.add_argument("a", type=float)
ap.add_argument("alpha", type=float)
args = ap.parse_args()


def check_grayimage(nameimage):
    """
    Check if image is monochrome (1 channel or 3 identical channels)
    """
    isGray = True
  
    im = PIL.Image.open(nameimage)
    if im.mode not in ("L", "RGB"):
        raise ValueError("Unsuported image mode for histogram computation")

    if im.mode == "RGB":
        pix = im.load()
        size = im.size
        for y in range(0, size[1]):
            if not isGray:
                break
            for x in range(0, size[0]):
                if not isGray:
                    break
                if (pix[x, y][0] != pix[x, y][1]) or \
                   (pix[x, y][0] != pix[x, y][2]):
                    isGray = False
    return isGray

def histogram(im, value_max=None, size=(256, 128), margin=10, padding=3):
    """
    Return an image object displaying the histogram of the input.
    Image mode supported : L, RGB
    """
    if im.mode not in ('L', 'RGB'):
        raise ValueError("Unsuported image mode for histogram computation (mode = {0})".format(im.mode))

    width_border = padding * 2 + 2 + size[0]
    height_border = padding * 2 + 3 + size[1]

    if im.mode == 'L':
        data = im.histogram()
        if not value_max:
            value_max = max(data)
        out = Image.new('L', (width_border, height_border), 255)
        draw = ImageDraw.Draw(out)
        draw_histo(draw, data, 'L', value_max, xy=(padding + 1, padding + 1), size=size, padding=padding)
        del draw
        out.value_max = value_max
        return out

    # polycolor in one frame
    if im.mode == 'RGB':
        data = im.histogram()
        # append gray intensity to data
        im_grey = im.convert("L")
        data += im_grey.histogram()
        if not value_max:
            value_max = max(data)
        out = Image.new('RGB', (width_border, height_border * 4 + margin * 3), (255, 255, 255))
        draw = ImageDraw.Draw(out)
        y = padding + 1
        for channel in ('R', 'G', 'B', 'I'):
            draw_histo(draw, data, channel, value_max, xy=(padding + 1, y), size=size, padding=padding)
            y += margin + height_border
        del draw
        del im
        out.value_max = value_max
        return out

def draw_histo(draw, data, channel, value_max, xy=(0, 0), size=(256, 128), padding=3):
    """
    Draws the histogram of a channel on an image.
    """
    dic = {
        # channel: (data_index, color, color_full, grid)
        'R': {'index': 0, 'color': (255, 0, 0), 'full': (0, 0, 0), 'grid': (192, 192, 192)},
        'G': {'index': 256, 'color': (0, 255, 0), 'full': (0, 0, 0), 'grid': (192, 192, 192)},
        'B': {'index': 512, 'color': (0, 0, 255), 'full': (0, 0, 0), 'grid': (192, 192, 192)},
        'I': {'index': 768, 'color': (128, 128, 128), 'full': (0, 0, 0), 'grid': (192, 192, 192)},
        'L': {'index': 0, 'color': 128, 'full': 0, 'grid': 192}
    }
    for y in (xy[1] + size[1] / 3, xy[1] + 2 * size[1] / 3): # horizontal grid
        draw.line((xy[0], y, xy[0] + size[0] - 1, y), fill=dic[channel]['grid'])
    for x in (xy[0] + size[0] / 3, xy[0] + 2 * size[0] / 3): # vertical grid
        draw.line((x, xy[1], x, xy[1] + size[1]), fill=dic[channel]['grid'])
    draw.rectangle(
        [(xy[0] - padding - 1, xy[1] - padding - 1), (xy[0] + size[0] + padding, xy[1] + size[1] + padding + 1)],
        outline=dic[channel]['full']
    ) # border

    sy = size[1] / float(value_max) # step y
    # loop on all pixelx of width
    for dx in range(0, size[0]):
        # index in data
        i = int(dx*256/float(size[0]))
        value = data[dic[channel]['index']+i]
        x = xy[0] + dx
        if not value:
            continue
        if value > value_max:
            color = dic[channel]['full']
            y2 = xy[1]
        else:
            color = dic[channel]['color']
            y2 = int(round(xy[1] + size[1] - (value * sy)))
        if y2 != xy[1] + size[1]:
            draw.line((x, xy[1] + size[1], x, y2), color)


p1 = ['poisson_lca', '-t', str(args.T), '-a', str(args.a), '-b', str(args.alpha), 'input_0.png', 'output_normI.png', 'output_darkI.png', 
                           'output_powerI.png', 'output_normRGB.png', 'output_darkRGB.png', 'output_powerRGB.png']
subprocess.run(p1)


#check if input image is monochrome
isGray = check_grayimage('input_0.png')

#Compute histograms of images
im0 = PIL.Image.open('input_0.png')
imI = PIL.Image.open('output_normI.png')
imdI = PIL.Image.open('output_darkI.png')
impI = PIL.Image.open('output_powerI.png')
if not isGray:
    with open('algo_info.txt', 'w') as file:
        file.write("notisGray=1")
    imRGB = PIL.Image.open('output_normRGB.png')
    imdRGB = PIL.Image.open('output_darkRGB.png')
    impRGB = PIL.Image.open('output_powerRGB.png')


#compute maximum of histogram values
dataH = im0.histogram()
value_max = max(dataH)

#draw all the histograms using the same reference maximum
im0 = histogram(im0, value_max=value_max)
im0.save('input_0_hist.png')
imI = histogram(imI, value_max=value_max)
imI.save('output_normI_hist.png' )
imdI = histogram(imdI, value_max=value_max)
imdI.save('output_darkI_hist.png' )
impI = histogram(impI, value_max=value_max)
impI.save('output_powerI_hist.png' )
if not isGray:
    impRGB = histogram(imRGB, value_max=value_max)
    impRGB.save('output_normRGB_hist.png' )
    imdRGB = histogram(imdRGB, value_max=value_max)
    imdRGB.save('output_darkRGB_hist.png' )
    impRGB = histogram(impRGB, value_max=value_max)
    impRGB.save('output_powerRGB_hist.png' )