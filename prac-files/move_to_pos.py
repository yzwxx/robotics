import math

def compute_coord(x, y, wx, wy):
	dx = wx - x
	dy = wy - y
	alpha = math.atan2(dy, dx)
	dist = math.sqrt(dx*dx + dy*dy)
	return (alpha, dist)

def compute_angle_turn(curr_angle, dest_angle):
	angle_diff = dest_angle - curr_angle
	if angle_diff > math.pi:
		angle_diff = -(math.pi*2 - angle_diff)
	if angle_diff < -math.pi:
		angle_diff = math.pi*2 + angle_diff
	return angle_diff

def debug_move_to_pos(curr_x,curr_y, curr_alpha, alpha, dist):
	new_x = curr_x + dist*math.cos(curr_alpha + alpha)
	new_y = curr_y + dist*math.sin(curr_alpha + alpha)
	new_alpha = curr_alpha + alpha
	print "New x, new y, new angle:"
	print new_x
	print new_y
	print new_alpha 
	
def debug(x,y,alpha,new_x, new_y):
	(new_alpha, dist) = compute_coord(x,y,new_x,new_y)
	angle_diff = compute_angle_turn(alpha, new_alpha)
	debug_move_to_pos(x,y,alpha, angle_diff, dist)
	
debug(6,9,0,-1,-1)
