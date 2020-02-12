# About

This repository contains the code and report for a CS course at Ecole Centrale de Lyon.
It mainly consists in a raytracer code written from scratch in C++.

# Results

The work is incremental. Each step consists in a defined set of features working.
The code for a given step correspond to a git tag.

## Step 1: Geometry detection

![step 1 render: white disc on black background](./result/step1.png)

## Step 2: Light diffusion

![step 2 render: white ball on balck background](./result/step2.png)

## Step 3: Multiple object handling and color handling

![step 3 render: white ball in a "room" made of red, blue and green balls](./result/step3.png)

"Step" 3.5: refactor code about material and intersection. Add a basic invert gamma function to reduce constrast.
![step 3.5 render: same as 3.5 with less constrast](./result/step3.5.png)

## Step 4: Add shadows

![step 4 render: previous scene with an additional litle ball projecting its shadow on the main one](./result/step4.png)

## Step 5: Add mirror like surfaces

![step 5 render: previous scene but the main ball is a spheric mirror](./result/step5.png)

## Step 6: Add transparent material

![step 6 render: previous scene but the main ball is made of glass](./result/step6.png)

"Step" 6.5: added stocastic raytracing for partial reflexions at transparent interfaces.
![step 6.5 render: two big balls side by side, made of glass and water, one little mirror ball](./result/step6.5.png)

## Step 7: Add indirect lighting and extended light sources

![step 7 render: left ball s a mirror, right is diffuse, shadows are smooth and not fully black](./result/step7.png)

"Step" 7.5: added depht-of-field effect. A lot of noise appear, it decrease with more ray launched (seed [high resolution render](./result/highres_3.png)).
![step 7.5 render: right ball is blurry and very close, background is blurry](./result/step7.5.png)
