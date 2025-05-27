import glob

from   PIL import Image

files = glob.glob('sers/ser_*.jpg')

for _file in files:
  img   = Image.open(_file)
  img   = img.convert("RGBA")

  datas = img.getdata()

  newData = []

  for item in datas:
    if (item[0] < 30) and (item[1] < 30) and (item[2] < 30):
        newData.append((0, 0, 0, 0))

    else:
        newData.append(item)

  img.putdata(newData)

  bits = _file.split('.')
  
  img.save(bits[0] + 'a.png', "PNG")
