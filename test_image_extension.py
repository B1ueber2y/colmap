import pycolmap

image = pycolmap.create_test_image()
print(type(image.extension))
print(image.extension.test_field)

