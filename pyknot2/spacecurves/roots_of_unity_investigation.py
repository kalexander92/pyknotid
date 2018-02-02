vk_list = ['2.1']

for i in range(1,8):
    knot_name = '3.' + str(i)
    vk_list.append(knot_name)

for i in range(1,109):
    knot_name = '4.' + str(i)
    vk_list.append(knot_name)


def root(number):
    root_unity = n.exp( 2 * n.pi * 1.j / number)
    root_unity = root_unity.real.round(3) + root_unity.imag.round(3)*1.j
    return root_unity


x, y, s, t = sym.symbols('x y s t')

#first_root = (-1)
#second_root = (-0.5 + (n.sqrt(3)/2)*1.j)
#third_root = (1.j)

#roots = [first_root, second_root, third_root]


#calc_poly = vwp.web_alex_from_gauss(knot).subs(y, x/y)
#web_poly = vwp.web_alex_poly_xy(knot)

roots_comparison = []

roots = [[2,3], [2,4], [3,4]]

for knots in vk_list:
    knots_info = [knots]
    #calc_poly = vwp.web_alex_from_gauss(knots).subs(y, x/y)
    calc_poly = vwp.web_alex_from_gauss(knots).subs([[s,x],[t,y]])
    web_poly = vwp.web_alex_poly(knots).subs([[s,x],[t,y]])
    for entries in roots:
        i, j = entries
        absolute_calc_poly = n.abs(calc_poly.subs([[x,root(i)],[y,root(j)]]).expand())
        absolute_web_poly = n.abs(web_poly.subs([[x,root(i)],[y,root(j)]]).expand())
        print 'i =', i, 'j =', j, root(i), root(j)
        print absolute_calc_poly - absolute_web_poly
        print absolute_calc_poly
        print absolute_web_poly
        print '\n'
        knots_info.append(absolute_calc_poly)
        roots_comparison.append([knots, i, j, abs(absolute_calc_poly - absolute_web_poly) < 0.01, absolute_calc_poly, absolute_web_poly])
    #roots_comparison.append(knots_info)


#for knots in vk_list:
#    calc_poly = vwp.web_alex_from_gauss(knots).subs(y, x/y)
#    web_poly = vwp.web_alex_poly_xy(knots)
#    for entries in roots:
#        i, j = entries
##for i in range(1,5):
##    for j in range(1,5):
#        absolute_calc_poly = n.abs(calc_poly.subs([[x,root(i)],[y,root(j)]]).expand())
#        absolute_web_poly = n.abs(web_poly.subs([[x,root(i)],[y,root(j)]]).expand())
#        print 'i =', i, 'j =', j, root(i), root(j)
#        print absolute_calc_poly - absolute_web_poly
#        print absolute_calc_poly
#        print absolute_web_poly
#        print '\n'
#        roots_comparison.append([knots, i, j, abs(absolute_calc_poly - absolute_web_poly) < 0.01, absolute_calc_poly, absolute_web_poly])


#for i in range(1,6):
#    for j in range(1,6):
#        for k in range(1,6):
#            for l in range(1,6):
#                print 'i =', i, 'j =', j, 'k = ', k, 'l = ', l
#                print abs((calc_poly.subs([[x,root(i,k)],[y,root(j,l)]]) - web_poly.subs([[s,root(i,k)],[t,root(j,l)]])).expand()) < 0.01
#                #print abs(calc_poly.subs([[x,root(i,k)],[y,root(j,l)]]).expand())
#                print abs(web_poly.subs([[s,root(i,k)],[t,root(j,l)]]).expand())
#                print '\n'


#for i in range(3):
#    for j in range(3):
#        print 'i =', roots[i], 'j =', roots[j]
#        print (calc_poly.subs([[x,roots[i]],[y,roots[j]]]) - web_poly.subs([[s,roots[i]],[t,roots[j]]])).expand()
#        #print calc_poly.subs([[x,roots[i]],[y,roots[j]]])
#        #print web_poly.subs([[x,roots[i]],[y,roots[j]]])
