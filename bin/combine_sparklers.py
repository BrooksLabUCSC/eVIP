import os
from PIL import Image, ImageDraw, ImageFont

def run(dir,out_image):
    fnt = ImageFont.truetype('/Library/Fonts/Arial.ttf', 20)


    #get mut folders in the eVIPP output directories
    mut_folders = [dI for dI in os.listdir(dir) if os.path.isdir(os.path.join(dir,dI))]

    sparkler_file_list = []
    for i in mut_folders:
        wtspec_sparkler = dir+"/"+i+"/wt_specific/eVIPP_sparkler_plots/"+i+"_spark_plots.png"
        mutspec_sparkler = dir+"/"+i+"/mutation_specific/eVIPP_sparkler_plots/"+i+"_spark_plots.png"
        sparkler_file_list.append([wtspec_sparkler,mutspec_sparkler])

    blank_image = Image.new('RGB', (2*640 +10, len(sparkler_file_list)*480 + 100),color = 'white')

    y_count = 0

    for l in sparkler_file_list:

        mut = l[0].split("/")[-4]

        d = ImageDraw.Draw(blank_image)
        d.text((10,y_count), mut, font=fnt, fill=(0, 0, 0))

        #WT-specific
        try:
            image_wt = Image.open(l[0])
            blank_image.paste(image_wt, (0,20+y_count))
            d = ImageDraw.Draw(blank_image)
            d.text((10,30+y_count), "WT-specific", font=fnt, fill=(0, 0, 0))

        except:
            pass

        #Mutation-specific
        try:
            image_mut = Image.open(l[1])
            blank_image.paste(image_mut, (650,20+y_count))
            d = ImageDraw.Draw(blank_image)
            d.text((650,30+y_count), "Mutation-specific", font=fnt, fill=(0, 0, 0))

        except:
            pass


        y_count +=530

    blank_image.save(out_image)

if __name__ == "__main__": main()
