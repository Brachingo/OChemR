from PIL import Image

# Load the images

white = Image.open("symbols/white.png")
plus_sign = Image.open("symbols/plus/plus1.png")

#Resize the plus sign if needed
#plus = plus.resize((white.width, white.height))
#plus_size = (50, 50)
#plus_sign = plus_sign.resize(plus_size)

# Calculate the position to center the plus sign
x_pos = (white.width - plus_sign.width) // 2
y_pos = (white.height - plus_sign.height) // 2

# Paste the plus sign on the white image

white.paste(plus_sign, (x_pos, y_pos), plus_sign)

# Save the image
white.save("output/plus.png")

from PIL import Image, ImageDraw

# Create a white background
img = Image.new('RGB', (200, 200), color='white')
draw = ImageDraw.Draw(img)

# Draw a plus sign
plus_color = 'black'
thickness = 10
center = img.width // 2

# Horizontal line
draw.line([(center - 50, center), (center + 50, center)], 
          fill=plus_color, width=thickness)

# Vertical line
draw.line([(center, center - 50), (center, center + 50)], 
          fill=plus_color, width=thickness)

# Save/show
img.save('created_plus.png')
img.show()