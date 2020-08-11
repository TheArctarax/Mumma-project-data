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
    draw.text((1200, 280), "d = {} Mpc".format(str(image)),(0,0,0),font=font)
    img.save("./processed_{}.png".format(str(image)))

