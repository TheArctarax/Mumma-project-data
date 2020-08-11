from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
import numpy as np

start_val = 110
end_val = 490
val_spacing = 10

im_array = np.arange(start_val, end_val + val_spacing, val_spacing)

for image in im_array:
    img = Image.open("./distance_{}_corner.png".format(str(image)))
    draw = ImageDraw.Draw(img)
    # font = ImageFont.truetype(<font-file>, <font-size>)
    font = ImageFont.truetype("./times_new_roman.ttf",180)
    # draw.text((x, y),"Sample Text",(r,g,b))
    x, y = 1150, 280
    draw.text((x, y), "d  = {} Mpc".format(str(image)),(0,0,0),font=font)
    # subscripts
    w = draw.textsize("d", font=font)[0]
    x += w
    y += font.getoffset("d")[1]

    font = ImageFont.truetype("./times_new_roman.ttf",100)
    h2 = draw.textsize("L", font=font)[1]
    y += h2 / 1.8

    draw.text((x, y), "L",(0,0,0),font=font)

    img.save("./processed_{}.png".format(str(image)))

