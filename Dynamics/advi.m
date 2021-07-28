function ad=advi(v,i)

ad=[skewop(v(1:3,i)) zeros(3,3);skewop(v(4:6,i)) skewop(v(1:3,i))];

end