function sita = GetAngle(p1, p2)
    temp = sum(p1.*p2)/(sqrt(sum(p1.^2))*sqrt(sum(p2.^2)));
    sita = acos(temp);
end