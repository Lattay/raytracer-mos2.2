;; monte_carlo.scm
;;
;; Scheme Chicken v5
;;
;; Compute \int_{-\pi/2}^{pi/2} \cos^{30}(x) dx
;; as an exercise
;; Results are clearly off but the code seems ok.
;; It may be a problem with the RNG being biased

(import chicken.random)

(define pi 3.141592654)

;; Produce a gaussian distributed number
(define (box-muller sigma mu)
  (let ((u1 (pseudo-random-real))
        (u2 (pseudo-random-real)))
    (+ mu
       (* sigma
          (sqrt (* -2 (log u1)))
          (cos (* 2 pi u2))))))

(define (box-muller2 sigma mu)
  (let ((u1 (pseudo-random-real))
        (u2 (pseudo-random-real)))
    (values 
      (+ mu
         (* sigma
            (sqrt (* -2 (log u1)))
            (cos (* 2 pi u2))))
      (+ mu
         (* sigma
            (sqrt (* -2 (log u1)))
            (sin (* 2 pi u2)))))))
 
;; Gaussian distribution of mean 1 and deviation 1
(define (gauss sigma x)
  (/ (exp (/ (* x x) (* -2 sigma sigma))) (* sigma (sqrt (* 2 pi)))))

;; compute \int_{-\pi/2}^{pi/2} \cos^{30}(x) dx with n samples
(define (monte-carlo-1 sigma n)
  (let sum ((s 0) (k 0))
    (if (= k n)
        (/ s n)
        (sum
          (+ s
             (let ((x (box-muller sigma 0)))
               (/ (expt (cos x) 30) (gauss sigma x))))
          (add1 k)))))

;; compute \int_{-\pi/2}^{pi/2} \int_{-\pi/2}^{pi/2} \int_{-\pi/2}^{pi/2} \cos^{30}(xyz) dxdydz with n samples
(define (monte-carlo-2 sigma n)
  (let sum ((s 0) (k 0))
    (if (= k n)
        (/ s n)
        (sum
          (+ s
             (let-values (((x y)  (box-muller2 sigma 0))
                          ((z _)  (box-muller2 sigma 0)))
               (/ (expt (cos (* x y z)) 30)
                  (gauss sigma x) (gauss sigma y) (gauss sigma z))))
          (add1 k)))))

(display (monte-carlo-1 0.25 1000000))
(newline)
(display (monte-carlo-2 0.25 1000000))
(newline)
