from PIL import Image, ImageDraw, ImageFont

# Load image
img = Image.open("MPM_Step1.png").convert("RGB")

# Create drawing context
draw = ImageDraw.Draw(img)

# Load font (you can use any .ttf file or default)
#font = ImageFont.truetype("arial.ttf", size=24)  # Or use 
font = ImageFont.load_default()

# Add text
draw.text((10, 10), "(a) Case A", font=font, fill=(0, 0, 0))  # Position, text, font, color

# Save or use this image later
img.save("labeled_image1.png")



