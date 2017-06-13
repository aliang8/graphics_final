import math

e = 120
V = [0, 0, 1]

def calculate_normal(polygons, i):

    A = [0, 0, 0]
    B = [0, 0, 0]
    N = [0, 0, 0]

    A[0] = polygons[i+1][0] - polygons[i][0]
    A[1] = polygons[i+1][1] - polygons[i][1]
    A[2] = polygons[i+1][2] - polygons[i][2]

    B[0] = polygons[i+2][0] - polygons[i][0]
    B[1] = polygons[i+2][1] - polygons[i][1]
    B[2] = polygons[i+2][2] - polygons[i][2]

    N[0] = A[1] * B[2] - A[2] * B[1]
    N[1] = A[2] * B[0] - A[0] * B[2]
    N[2] = A[0] * B[1] - A[1] * B[0]

    return N

def sub_vect ( v1, v2 ):
    return [v1[i] - v2[i] for i in range(3)]

def scale ( v, s ):
    return [s * i for i in v]

def magnitude ( v ):
    return math.sqrt(sum(component**2 for component in v))

def dot_prod ( v1, v2 ):
    v1 = normalize(v1)
    v2 = normalize(v2)
    p = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    return max(0, p)

def normalize ( v ):
    mag = magnitude(v)
    return [component / mag for component in v ]

# I_ambient ( i_a = A * k_a )
def i_a ( L_c, k_a ):
    return L_c * k_a

# I_diffuse ( i_d = L * K_d * max(0, N * L) ) <- Lambert's Law
# where L is the point light source
def i_d ( L_c, k_d, N, L ):
    return L_c * k_d * max(0, dot_prod(N, L))

# I_specular ( i_s = L * k_s * max(0, R * V) ^e )
e = 120
def i_s ( L_c, k_s, N, L ):
    p = scale(scale(N, 2), dot_prod(N, L))
    R = sub_vect(p, L)

    return L_c * k_s * (max(0, dot_prod(R, V)) ** e)

def intensity( i_a, i_d, i_s ):
    sum = i_a + i_d + i_s
    return min(255, int(round(sum)))

def calc_lighting( N, constants, source ):
    L = source['location']

    a_r, d_r, s_r = constants['red']
    a_b, d_b, s_b = constants['blue']
    a_g, d_g, s_g = constants['green']

    L_r, L_g, L_b = source['color']


    red = intensity(i_a(L_r, a_r), i_d(L_r, d_r, N, L), i_s(L_r, s_r, N, L))
    green = intensity(i_a(L_g, a_g), i_d(L_g, d_g, N, L), i_s(L_g, s_g, N, L))
    blue = intensity(i_a(L_b, a_b), i_d(L_b, d_b, N, L), i_s(L_b, s_b, N, L))

    return [red, green, blue]


def total_lighting ( N, constants, source ):
    color = [0, 0 ,0]

    for s in source:
        light = calc_lighting(N, constants, s)
        color[0] += light[0]
        color[1] += light[1]
        color[2] += light[2]

    return [min(255, i) for i in color]

def calc_average(normals):
    x = y = z = 0.0

    for n in normals:
        x += n[0]
        y += n[1]
        z += n[2]

    length = len(normals)
    return [ x/length, y/length, z/length ]

# get vertex normals for phong shading
def get_vertex_normals(matrix):
    v_n = {}

    point = 0
    while point < len(matrix) - 2:
        normal = calculate_normal(matrix, point)[:]

        v_1 = (int(matrix[point][0]), int(matrix[point][1]), matrix[point][2])
        v_2 = (int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2])
        v_3 = (int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2])
        
        if v_1 in v_n:
            v_n[v_1].append(normal)
        else:
            v_n[v_1] = [ normal ]

        if v_2 in v_n:
            v_n[v_2].append(normal)
        else:
            v_n[v_2] = [ normal ]

        if v_3 in v_n:
            v_n[v_3].append(normal)
        else:
            v_n[v_3] = [ normal ]

        point += 3

    for key in v_n:
        v_n[key] = calc_average(v_n[key])

    return v_n
