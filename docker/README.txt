## Build image from your Dockerfile instruction
docker build -t organmatch -f docker/Dockerfile .

## Run the app and test if it works well
docker run -p 3839:3839 organmatch:latest

## You can also interactively run the image and test all the functions 
docker run -it --user=root organmatch:latest bash

## Instide the image, try `R` and loading all libraries
R --vanilla + ENTER
library(...) # just like you would in R 
ctrl + d
exit # to exit interatcive shell
# This is also the best way to diagnose why installation of some packages or libraries failed 
# and test potential solutions

## Save image to tar file for porting
docker save -o docker/organmatch.tar organmatch:latest

## In the new host.. load the tarred image
docker load -i docker/organmatch.tar
