import pycolmap

image = pycolmap.create_test_image()
extension = image.get_extension()
print(type(extension))
print(extension.test_field)

