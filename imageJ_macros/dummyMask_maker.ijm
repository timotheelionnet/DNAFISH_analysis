//Run this macro to create a dummy mask 
//that contains your whole FOV

//Set params of the images you are working with
imXYDim = 1200; 
depth = "16-bit";

newImage("dummyMask", "Black", imXYDim, imXYDim, depth)
setColor("white");
drawRect(1, 1, 1198, 1198);
fill();


