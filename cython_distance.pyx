cpdef double calculateDistance((double,double) point, (double,double) ppoint, (double,double) immuneCell):
    cdef double distance
    cdef double normalx, normaly, grad, dbpoints, dimmuneToPoint, dimmuneToPoint1, dimmuneToPoint2, difference
    if point[0]-ppoint[0]==0.0 and point[1]-ppoint[1]==0.0:
        distance = ((immuneCell[0]-point[0])**2+(immuneCell[1]-point[1])**2)**0.5
    else:
        if point[0]-ppoint[0]==0.0:
            normaly=0.0
            normalx=1.0
            grad=0.0
        elif point[1]-ppoint[1]==0.0:
            normaly=1.0
            normalx=0.0
            grad=0.0
        else:
            grad = (point[1]-ppoint[1])/(point[0]-ppoint[0])
            normalx = 1.0/(1.0**2+1/grad**2)**0.5
            normaly = -(1.0/grad)/(1.0**2+1.0/grad**2)**0.5

        dbpoints = ((ppoint[0]-point[0])**2+(ppoint[1]-point[1])**2)**0.5
        dimmuneToPoint1 = ((immuneCell[0]-point[0])**2+(immuneCell[1]-point[1])**2)**0.5
        dimmuneToPoint2 = ((immuneCell[0]-ppoint[0])**2+(immuneCell[1]-ppoint[1])**2)**0.5
        dimmuneToPoint=max(dimmuneToPoint1,dimmuneToPoint2)

        distance = abs(normalx*(point[0]-immuneCell[0]) + normaly*(point[1]-immuneCell[1]))
        if dbpoints**2+distance**2 < dimmuneToPoint**2 :
            distance = min(dimmuneToPoint1,dimmuneToPoint2)
        if distance>dimmuneToPoint:
            difference = distance-dimmuneToPoint
            if difference>1e-7:
                print("Differnce:      "+ str(difference))
                print("point:          "+str(point))
                print("previous point: "+str(ppoint))
                print("immunceCell:    "+str(immuneCell))
                print("dbpoints:       "+str(dbpoints))
                print("dimmuneToPoint: "+str(dimmuneToPoint))
                print("distance:       "+str(distance))
                print("normalx:        "+str(normalx))
                print("normaly:        "+str(normaly))
                print("grad:           "+str(grad))
        
    return distance

cpdef double calculateSDistance((double,double) point,(double,double) immuneCell):
    cdef double distance						  
    distance = ((immuneCell[0]-point[0])**2+(immuneCell[1]-point[1])**2)**0.5
    return distance

if __name__=="__main__":
    point=(14905.86, 7369.56)
    ppoint=(14899.83, 7375.59)
    immuneCell=(12346.3, 4810.0)

    distance = calculateDistance(point,ppoint,immuneCell)
