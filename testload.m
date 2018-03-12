  


cac = data();




function    cac = data()
        str = fileread( 'testdata1.txt' );
        str = strrep( str, ' ', '' );
        cac = textscan( str, '%f%f%f%f', 'Delimiter', ',' ); 
end