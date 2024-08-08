# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:53:38 2024

@author: dpqb1
"""

from PIL import Image
import os

def remove_whitespace(input_dir, output_dir,threshold=100):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".png"):
            # Open the image
            image_path = os.path.join(input_dir, filename)
            img = Image.open(image_path)
            
            # Convert the image to grayscale
            img_gray = img.convert("L")

            # Apply a threshold to the image
            img_thresh = img_gray.point(lambda p: p < threshold and 255)

            # Find bounding box of non-white pixels
            bbox = img_thresh.getbbox()
            
            # Remove whitespace
            img = img.crop(bbox)

            # Save the image
            output_path = os.path.join(output_dir, filename)
            img.save(output_path)

            print(f"Whitespace removed from {filename}")

# Replace these paths with your input and output directories
input_directory = "C:\\Users\\dpqb1\\Documents\\ONR"
output_directory = "C:\\Users\\dpqb1\\Documents\\ONR"

remove_whitespace(input_directory, output_directory)
