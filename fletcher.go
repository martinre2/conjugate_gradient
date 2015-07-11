package main

import (
	"fmt"
	"math"
)

const delta_x float64 = 0.00001
const delta float64 = 0.00001

var K uint = 0

func function_test(vars []float64) float64 {
	x1 := vars[0]
	x2 := vars[1]

	return x1 - x2 + (2 * math.Pow(x1, 2)) + (2 * (x1 * x2)) + math.Pow(x2, 2)
}
func function_test2(vars []float64) float64 {
	x1 := vars[0]
	x2 := vars[1]
	x3 := vars[2]
	x4 := vars[3]

	m1 := func() float64 { return 10 * (x2 - math.Pow(x1, 2)) }
	m2 := func() float64 { return math.Pow(1-x1, 2) }
	m3 := func() float64 { return 90 * (x4 - math.Pow(x3, 2)) }
	m4 := func() float64 { return math.Pow(1-x3, 2) }
	m5 := func() float64 { return 10 * math.Pow(x2+x4-2, 2) }
	m6 := func() float64 { return 0.1 * (x2 - x4) }

	return math.Pow(m1(), 2) + m2() + m3() + m4() + m5() + m6()
}

func calc_gradient(vars []float64) []float64 {
	s_direction := make([]float64, len(vars))

	for i := range vars {
		s_direction[i] = -1.0 * gradient(vars, uint(i))
	}
	return s_direction
}

func gradient(vars []float64, i uint) float64 {
	var m1, m2 float64

	b := make([]float64, len(vars))
	copy(b, vars)

	b[i] = b[i] + delta_x
	m1 = function_test(b)

	b[i] = b[i] - (2 * delta_x)
	m2 = function_test(b)

	return (m1 - m2) / (2 * delta_x)
}

func calc_residual_vector(vars []float64) []float64 {
	residue := make([]float64, len(vars))

	if K == 0 {
		for i, _ := range vars {
			residue[i] = gradient(vars, uint(i))
		}

	} else {

	}
	return residue
}

func golden_section(a []float64, b []float64, tolerance float64) []float64 {
	L := minus(b, a)

	for {
		x1 := add(a, mul(L, 0.618))
		x2 := minus(b, mul(L, 0.618))

		if function_test(x1) < function_test(x2) {
			a = x2
		} else {
			b = x1
		}
		L = minus(b, a)

		if norm(L) < tolerance {
			break
		}
	}

	return div(add(a, b), 2)
}

func acotamiento(vars []float64, delta float64) ([]float64, []float64) {
	k := 0
	xs := make([]*[]float64, 3)
	var in_delta float64
	in_delta = delta

	b := make([]float64, len(vars))
	copy(b, vars)

	xs[k] = &b
	i := 0

	for {
		x := function_test(*xs[i])
		x_minus_delta := resolve_with_delta(*xs[i], (in_delta * -1.0))
		x_plus_delta := resolve_with_delta(*xs[i], (in_delta * 1.0))

		//fmt.Println(x, x_minus_delta, x_plus_delta)

		if x <= x_minus_delta && x >= x_plus_delta {
			//fmt.Println("Case 1")
			in_delta = (in_delta) * 1.0
		} else if (x >= x_minus_delta) && (x <= x_plus_delta) {
			//fmt.Println("Case 2")
			in_delta = in_delta * -1.0
		} else if (x <= x_minus_delta) && (x <= x_plus_delta) {
			//fmt.Println("Case 3")
			break
		}
		k++
		i = k

		x_next := apply_delta(b, math.Pow(2, float64(k-1))*in_delta)

		if i > 2 {
			i = 2
			remove(xs, 0)
			xs[i] = &x_next
		} else {
			xs[i] = &x_next
		}

		eval_x_next := function_test(x_next)

		//fmt.Println(x, eval_x_next)

		if eval_x_next < x && k > 2 {
			break
		}
		//fmt.Println(xs)

	}

	//fmt.Println(xs, k)

	return min(*xs[0], *xs[2]), max(*xs[0], *xs[2])
}

func resolve_with_delta(vars []float64, delta float64) float64 {

	fmt.Println("with delta dlt:", delta)
	fmt.Println("with delta:", vars)
	b := make([]float64, len(vars))
	copy(b, vars)

	for i := range b {
		b[i] += delta
	}
	fmt.Println("with delta:", b)

	return function_test(b)
}

func apply_delta(vars []float64, delta float64) []float64 {

	b := make([]float64, len(vars))
	copy(b, vars)

	for i := range b {
		b[i] += delta
	}

	return b
}

func remove(arr []*[]float64, item int) {
	arr = append(arr[:item], arr[item+1:]...)
}

func min(a []float64, b []float64) []float64 {
	r := make([]float64, len(a))

	for i := range a {
		if a[i] < b[i] {
			r[i] = a[i]
		} else {
			r[i] = b[i]
		}
	}

	return r
}

func max(a []float64, b []float64) []float64 {
	r := make([]float64, len(a))

	for i := range a {
		if a[i] > b[i] {
			r[i] = a[i]
		} else {
			r[i] = b[i]
		}
	}

	return r
}

func minus(a []float64, b []float64) []float64 {
	r := make([]float64, len(a))

	for i := range a {
		r[i] = a[i] - b[i]
	}

	return r
}

func mul(v []float64, scalar float64) []float64 {
	r := make([]float64, len(v))

	for i := range v {
		r[i] = v[i] * scalar
	}

	return r
}

func div(v []float64, scalar float64) []float64 {
	r := make([]float64, len(v))

	for i := range v {
		r[i] = v[i] / scalar
	}

	return r
}

func add(a []float64, b []float64) []float64 {
	r := make([]float64, len(a))

	for i := range a {
		r[i] = a[i] + b[i]
	}

	return r
}

func norm(v []float64) float64 {
	var sum float64
	for _, d := range v {
		sum += math.Pow(d, 2)
	}

	return math.Sqrt(sum)
}

func norm_gradient(vars []float64) float64 {

	var sum float64

	for i := range vars {
		sum += math.Pow(gradient(vars, uint(i)), 2)
	}
	return sum
}

func main() {
	/*vars := []float64{
		-3,
		-1,
		-3,
		-1,
	}*/

	vars := []float64{
		-1,
		1,
	}

	var k int = 0
	s_direction := calc_gradient(vars)
	fmt.Println("s0:", s_direction)

	a, b := acotamiento(vars, delta)
	lamda := function_test(golden_section(a, b, 0.00001))

	x_next := add(vars, mul(s_direction, lamda))
	x_prev := make([]float64, len(x_next))
	s_prev := make([]float64, len(s_direction))
	copy(x_prev, x_next)
	copy(s_prev, s_direction)

	for {
		k++
		fmt.Println("x", k, x_next, function_test(x_next))

		s_direction = add(calc_gradient(x_next), mul(s_prev, (norm_gradient(x_next)/norm_gradient(x_prev))))
		a, b = acotamiento(x_next, delta)
		lamda = function_test(golden_section(a, b, 0.00001))

		copy(x_prev, x_next)
		copy(s_prev, s_direction)

		x_next = add(vars, mul(s_direction, lamda))
	}
}
