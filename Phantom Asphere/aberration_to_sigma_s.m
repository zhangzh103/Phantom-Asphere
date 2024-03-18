function w=aberration_to_sigma(w,NA)
w(4)=w(4)/NA*2
w(1)=w(1)/NA*4
w(2)=w(2)/NA*3
w(3)=w(3)/NA*3+w(4)
end