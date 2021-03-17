good_d = []
with open("../data/good-d-50000.txt") as fin:
    for d in fin:
        good_d.append(int(d))

good_d = [abs(d) for d in good_d]
good_d.sort()

D_RANGE = 50000

plot_d = []
ct = 0
for d in good_d:
    if abs(d) < D_RANGE:
        ct += 1
        plot_d.append((d,ct))
    
    

plt = plot([])
plt += list_plot(plot_d, color='blue')

var('a','b')
model(x) = a*(x**(1/3))
sol = find_fit(plot_d, model)
show(sol)
g(x) = model(a = sol[0].rhs())
plt += plot(g(x),  x, [1,D_RANGE], color = 'blue')
f(x) = 3.0863*x^(1/3)
h(x) = 3.2932*x^(1/3)
plt += plot(h(x), x, [1,D_RANGE], color = 'purple', axes_labels = ['D', "Number of good $d$ with $|d| \leq D$"])
plt += plot(g(x), x, [1,D_RANGE], color = 'purple')

save(plt,"good-d.png")

