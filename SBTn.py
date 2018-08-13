from shapely.geometry import Polygon, Point
import math as m

def returnSoilType(point):
    Polygon1 = Polygon([(-1.0000000000000004, 0.011161060300725603), (-1.0000000000000004, 1.0156564873662957), (-0.7636363636363641, 1.004495427065567), (-0.5402597402597407, 0.9338087118276195), (-0.4025974025974033, 0.8742830568904006), (-0.2441558441558449, 0.7589521004495386), (-0.101298701298702, 0.6138583165400677), (0.04935064935064837, 0.42784064486125867), (0.1402597402597392, 0.24926368004960156), (0.21818181818181737, 0.05952565493721625), (0.22857142857142754, -2.606755244094819e-15), (-1.0000000000000004, 0.011161060300725603)])
    Polygon2 = Polygon([(0.23116883116883025, -2.606755244094819e-15), (0.22337662337662256, 0.04464424120291128), (0.3350649350649342, 0.09300883583940192), (0.5246753246753235, 0.21950085258099214), (0.6597402597402591, 0.3422725158890061), (0.7688311688311678, 0.4687645326305963), (0.8909090909090898, 0.6510618508758297), (0.9636363636363623, 0.8221981088203337), (0.9948051948051939, 0.8928848240582813), (0.9999999999999991, 0.011161060300725603), (0.23116883116883025, -2.606755244094819e-15)])
    Polygon3 = Polygon([(0.22337662337662256, 0.055805301503640006), (0.14545454545454461, 0.25670438691675407), (0.04935064935064837, 0.44272205859556313), (-0.0987012987012994, 0.6212990234072202), (-0.14285714285714346, 0.6733839714772866), (-0.03636363636363715, 0.7254689195473529), (0.10909090909090846, 0.8296388156874862), (0.2883116883116874, 0.9821733064641098), (0.45714285714285596, 1.164470624709343), (0.5948051948051936, 1.369090063556033), (0.6779220779220771, 1.5104634940319281), (0.7532467532467519, 1.6815997519764325), (0.7844155844155836, 1.803673849015651), (0.8415584415584407, 1.7671678809486848), (0.9246753246753237, 1.748566113780804), (0.9999999999999991, 1.7448457603472276), (0.9948051948051939, 0.90404588435901), (0.96103896103896, 0.82591846225391), (0.8805194805194794, 0.6547822043094059), (0.7636363636363623, 0.4762052394977488), (0.6545454545454539, 0.34599286932258233), (0.5220779220779213, 0.22322120601456838), (0.3324675324675317, 0.09672918927297817), (0.22337662337662256, 0.055805301503640006)])
    Polygon4 = Polygon([(-0.5168831168831174, 0.9300883583940432), (-0.3662337662337669, 0.9896140133312623), (-0.16623376623376684, 1.0975042629049716), (-0.02597402597402665, 1.183072391877224), (0.11948051948051863, 1.305844055185238), (0.27012987012986933, 1.4472174856611328), (0.4311688311688304, 1.648116571074247), (0.5142857142857133, 1.793210354983718), (0.5818181818181809, 2.0164315609982895), (0.6675324675324665, 1.893659897690275), (0.7272727272727262, 1.8415749496202083), (0.7840909090909076, 1.803673849015651), (0.751948051948051, 1.6815997519764327), (0.6766233766233756, 1.5123236707487162), (0.5928571428571421, 1.3681599751976392), (0.45454545454545414, 1.165400713067737), (0.2870129870129863, 0.9858936598976861), (0.10649350649350575, 0.8314989924042746), (-0.03896103896103975, 0.728259184622535), (-0.1525974025974035, 0.6733839714772866), (-0.24220779220779298, 0.7645326305999033), (-0.40129870129870215, 0.8770733219655826), (-0.5168831168831174, 0.9300883583940432)])
    Polygon5 = Polygon([(-0.9948051948051947, 1.4341962486436162), (-0.8964285714285722, 1.436521469539602), (-0.7795454545454552, 1.4532630599906948), (-0.5211038961038963, 1.514648891644702), (-0.35876623376623384, 1.5816152534490733), (-0.13019480519480497, 1.7155479770578164), (0.023051948051948434, 1.8420399937994067), (0.22564935064935154, 2.0671213765307663), (0.3568181818181828, 2.273600992094245), (0.4035714285714296, 2.366609827933649), (0.4490259740259752, 2.275461168811033), (0.5438311688311701, 2.100604557432952), (0.5811688311688323, 2.0178266935358815), (0.5136363636363648, 1.7932103549837193), (0.4306818181818195, 1.648581615253445), (0.2701298701298709, 1.4476825298403306), (-0.026623376623375883, 1.1833049139668226), (-0.16639610389610338, 1.0982018291737672), (-0.3660714285714284, 0.9907766237792547), (-0.5243506493506493, 0.9282281816772548), (-0.5405844155844157, 0.9352038443652102), (-0.764935064935065, 1.0072856921407487), (-0.9987012987012989, 1.016586575724689), (-0.9948051948051947, 1.4341962486436162)])
    Polygon6 = Polygon([(-0.9977272727272727, 2.197333746705933), (-0.9509740259740257, 2.201054100139509), (-0.8262987012987009, 2.2289567508913306), (-0.6977272727272721, 2.275461168811033), (-0.5496753246753241, 2.33684700046504), (-0.4224025974025968, 2.414044334211746), (-0.32045454545454477, 2.490311579600058), (-0.19772727272727186, 2.612153154549678), (-0.11201298701298612, 2.746085878158421), (-0.05357142857142749, 2.847465509223372), (0.0068181818181829446, 2.999999999999996), (0.14707792207792347, 2.9981398232832084), (0.16525974025974155, 2.91815222446132), (0.20292207792207928, 2.832584095489068), (0.2568181818181827, 2.6958611068051423), (0.31071428571428683, 2.5749496202139164), (0.3522727272727286, 2.4819407843745114), (0.40292207792207946, 2.3684700046504377), (0.35487012987013156, 2.272205859556654), (0.22435064935065108, 2.066656332351569), (0.021428571428572463, 1.8439001705161948), (-0.13149350649350555, 1.7169431095954073), (-0.3600649350649344, 1.5830103859866642), (-0.5220779220779217, 1.5160440241822928), (-0.7805194805194805, 1.4541931483490884), (-0.8977272727272728, 1.4383816462563899), (-0.9951298701298704, 1.4355913811812075), (-0.9977272727272727, 2.197333746705933)])
    Polygon7 = Polygon([(-0.9980519480519495, 2.9995349558207978), (0.0063311688311682435, 2.999999999999995), (-0.05405844155844208, 2.84769803131297), (-0.11266233766233846, 2.7463184002480183), (-0.19821428571428656, 2.612385676639276), (-0.3207792207792216, 2.491241667958451), (-0.423051948051949, 2.4142768563013437), (-0.5500000000000012, 2.3375445667338344), (-0.6980519480519491, 2.2766237792590243), (-0.8266233766233778, 2.2298868392497235), (-0.951136363636365, 2.202216710587501), (-0.9975649350649364, 2.1984963571539247), (-0.9980519480519495, 2.9995349558207978)])
    Polygon8 = Polygon([(0.14724025974026023, 2.9983723453728053), (0.6912337662337673, 2.9974422570144115), (0.6925324675324682, 2.885831654007126), (0.6866883116883125, 2.7463184002480183), (0.6769480519480531, 2.6151759417144578), (0.6639610389610402, 2.4877538366144734), (0.6483766233766244, 2.3407998759882136), (0.6282467532467546, 2.2170981243218053), (0.6024350649350663, 2.076422260114706), (0.5818181818181831, 2.018524259804676), (0.5441558441558454, 2.1013021237017466), (0.4495129870129879, 2.276391257169426), (0.40389610389610464, 2.369167570919232), (0.31120129870129953, 2.5756471864827106), (0.2576298701298707, 2.696558673073937), (0.20308441558441648, 2.8339792280266574), (0.16574675324675403, 2.919082312819713), (0.14724025974026023, 2.9983723453728053)])
    Polygon9 = Polygon([(0.6913961038961056, 2.9976747791040093), (0.9991883116883149, 2.9993024337311986), (0.9998376623376652, 1.744962021392026), (0.9247564935064962, 1.7489148969152006), (0.8414772727272757, 1.7677491861726802), (0.7274350649350683, 1.8430863432025981), (0.66688311688312, 1.896566423810256), (0.5824675324675352, 2.017942954580679), (0.6030844155844188, 2.076305999069905), (0.6288961038961067, 2.217446907456202), (0.649025974025977, 2.3409161370330116), (0.6642857142857168, 2.487405053480074), (0.6775974025974056, 2.615292202759256), (0.6873376623376657, 2.746667183382415), (0.6930194805194838, 2.8854828708727265), (0.6913961038961056, 2.9976747791040093)])

    DomainList = [Polygon1, Polygon2, Polygon3, Polygon4, Polygon5, Polygon6, Polygon7, Polygon8, Polygon9] 
    if point[0] == 0:
        point[0] = 0.000000001


    if point[1] == 0:
        point[1] = 0.000000001
    x = m.log10(point[0])
    y = m.log10(point[1])
    if x < -1:
        x = -1
    if x > 1:
        x = 1
    if y < 0:
        y = 1
    if y > 3:
        y = 3

    thisLocation = Point(x,y)
    count = 0
    Domain = 0
    for thisPoly in DomainList:
        count += 1
        if thisLocation.within(thisPoly):
            Domain = count

    return Domain

pa = 0.1013
pointx = 0.006
pointy = 0.204/pa
point = [pointx, pointy]
Domain = returnSoilType(point)
print(Domain)