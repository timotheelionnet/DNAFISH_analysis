function movieInfo2D = movie_info3D_to_2D(movieInfo3D)

for i=1:size(movieInfo3D,2)
    movieInfo2D(i).xCoord = movieInfo3D(i).xCoord;
    movieInfo2D(i).yCoord = movieInfo3D(i).yCoord;
    movieInfo2D(i).amp = movieInfo3D(i).amp;
    
    
end