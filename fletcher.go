package main

import (
	"fmt"
	"math"
	"math/rand"
)

const delta_x float64 = 0.00000001
const delta float64 = 0.00001
const e float64 = 0.00001
const e2 float64 = 0.00001
const e3 float64 = 0.00001

var debug bool = true

var K uint = 0

func function_test(vars []float64) float64 {
	x1 := vars[0]
	x2 := vars[1]

	return 2*math.Pow(x1, 2) + (x1 * x2) - math.Pow(x2, 2)
}

func ddwfunction_test(vars []float64) float64 {
	x1 := vars[0]
	x2 := vars[1]

	return 2*x1 - math.Pow(x2, 2)
}

func d4function_test(vars []float64) float64 {
	x1 := vars[0]
	x2 := vars[1]

	return x1 - x2 + (2 * math.Pow(x1, 2)) + (2 * (x1 * x2)) + math.Pow(x2, 2)
}

func dfunction_test(vars []float64) float64 {
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

	b[i] = vars[i] + delta_x
	m1 = function_test(b)

	copy(b, vars)
	b[i] = vars[i] - (delta_x)
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

func golden_section(vars []float64, s []float64, a float64, b float64, tolerance float64) float64 {
	//fmt.Println("Gold")
	vc := make([]float64, len(vars))
	copy(vc, vars)

	L := b - a

	for {
		x1 := a + (L * 0.618)
		x2 := b - (L * 0.618)

		if function_test(add(vc, mul(s, x1))) < function_test(add(vc, mul(s, x2))) {
			a = x2
		} else {
			b = x1
		}
		L = (b - a)

		if math.Abs(L) < tolerance {
			break
		}
	}

	return (a + b) / 2
}

func acotamiento(vars []float64, s []float64, rx float64, delta float64) (float64, float64) {
	//fmt.Println("Acot")
	k := 0
	xs := make([]float64, 3)
	var in_delta float64
	in_delta = delta

	b := make([]float64, len(vars))
	s_rx := make([]float64, len(vars))
	s_rx_plus_delta := make([]float64, len(vars))
	s_rx_minus_delta := make([]float64, len(vars))
	copy(b, vars)

	i := 0

	for {
		//fmt.Println("delta:", in_delta)
		//fmt.Println("x0:", rx)
		//fmt.Println("x:", b)
		//fmt.Println("s:", s)

		s_rx = mul(s, rx)
		s_rx_plus_delta = mul(s, rx+(math.Abs(in_delta)))
		s_rx_minus_delta = mul(s, rx-(math.Abs(in_delta)))

		//fmt.Println("rx", s_rx)
		//fmt.Println("rx+d", s_rx_plus_delta)
		//fmt.Println("rx-d", s_rx_minus_delta)

		x := function_test(add(b, s_rx))
		x_plus_delta := function_test(add(b, s_rx_plus_delta))
		x_minus_delta := function_test(add(b, s_rx_minus_delta))

		//fmt.Println(x, x_minus_delta, x_plus_delta)

		if x <= x_minus_delta && x >= x_plus_delta {
			//fmt.Println("Case 1")
			in_delta = math.Abs(in_delta) * 1.0
		} else if (x >= x_minus_delta) && (x <= x_plus_delta) {
			//fmt.Println("Case 2")
			in_delta = math.Abs(in_delta) * -1.0
		} else if (x <= x_minus_delta) && (x <= x_plus_delta) {
			//fmt.Println("Case 3")
			break
		}
		k++
		i = k

		if i > 2 {
			i = 2
			remove(xs, 0)
			xs[i] = rx
		} else {
			xs[i] = rx
		}

		//rx += math.Pow(2, float64(k-1)) * in_delta
		rx += in_delta

		eval_prev_delta := function_test(add(b, mul(s, xs[i])))
		eval_in_delta := function_test(add(b, mul(s, rx)))

		//fmt.Println("eval:", eval_in_delta, eval_prev_delta)

		if eval_in_delta < eval_prev_delta {
			break
		}
		//fmt.Println(xs)

	}

	//fmt.Println(xs, k)

	if k < 2 {
		return rx - math.Abs(in_delta), rx + math.Abs(in_delta)
	} else {
		return min(xs[0], xs[2]), max(xs[0], xs[2])
	}
	return rx - math.Abs(in_delta), rx + math.Abs(in_delta)
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

func remove(arr []float64, item int) {
	arr = append(arr[:item], arr[item+1:]...)
}

func min_array(a []float64, b []float64) []float64 {
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

func min(a float64, b float64) float64 {
	if a < b {
		return a
	} else {
		return b
	}
}

func max_array(a []float64, b []float64) []float64 {
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

func max(a float64, b float64) float64 {
	if a > b {
		return a
	} else {
		return b
	}
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
	return math.Sqrt(sum)
}

func main() {
	/*
		vars := []float64{
			-3,
			-1,
			-3,
			-1,
		}
	*/
	vars := []float64{
		0.5,
		1,
	}

	var k int = 0

	fmt.Println("=> x", k, vars, function_test(vars))

	s_direction := calc_gradient(vars)
	//fmt.Println("s0:", s_direction)

	a, b := acotamiento(vars, s_direction, rand.Float64(), delta)
	//a, b := acotamiento(vars, s_direction, -0.2, delta)
	lamda := golden_section(vars, s_direction, a, b, e)
	//fmt.Println("lamda", lamda)

	x_next := add(vars, mul(s_direction, lamda))
	x_prev := make([]float64, len(x_next))
	s_prev := make([]float64, len(s_direction))
	copy(x_prev, vars)
	copy(s_prev, s_direction)
	for {
		k++
		fmt.Println("=> x", k, x_next, function_test(x_next))

		//fmt.Println("GK+1", calc_gradient(x_next))
		//fmt.Println("NORM x_n", norm_gradient(x_next))
		//fmt.Println("NORM x_p", norm_gradient(x_prev))
		s_direction = add(calc_gradient(x_next), mul(s_prev, (norm_gradient(x_next)/norm_gradient(x_prev))))
		//fmt.Println("s", k, s_direction)
		a, b = acotamiento(x_next, s_direction, rand.Float64(), delta)
		lamda = golden_section(x_next, s_direction, a, b, e)
		//fmt.Println("lamda", lamda)

		copy(x_prev, x_next)
		copy(s_prev, s_direction)

		x_next = add(x_next, mul(s_direction, lamda))

		if norm(minus(x_next, x_prev))/norm(x_prev) < e2 {
			break
		} else if norm_gradient(x_next) < e3 {
			break
		}

	}
}
