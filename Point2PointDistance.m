function d = Point2PointDistance(Lon, Lat)
Re= 6378;
Lat1 = Lat(1);
Lat2 = Lat(2);
Lon1 = Lon(1);
Lon2 = Lon(2);

d = Re * acos((sin(Lat1) * sin(Lat2) + cos(Lat1) * cos(Lat2) * cos(Lon2 - Lon1)));

end