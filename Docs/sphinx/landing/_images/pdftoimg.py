from pdf2image import convert_from_path
import sys

filename_pdf=sys.argv[1]
filename_png=filename_pdf[0:-3]+"png"
print("The pdf file name is "+filename_pdf)
print("The png file name is "+filename_png)

images = convert_from_path(filename_pdf, dpi=300)
images[0].save(filename_png, "PNG")
