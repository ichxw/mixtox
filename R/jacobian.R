jacobian <- function(eq, x, paraHat){
	## Jacobian matrix calculation
	n <- length(x)
	mpara <- length(paraHat)
	Alpha <- paraHat[1]
	Beta <- paraHat[2]
	if (mpara == 3) Gamma <- paraHat[3]
	if (mpara == 4) {Gamma <- paraHat[3]; Delta <- paraHat[4]}
	if (mpara == 5) {Gamma <- paraHat[3]; Delta <- paraHat[4]; Epsilon <- paraHat[5]}
	
	jac <- matrix(rep(0, n * mpara), n, mpara)
	
	jacFun <- switch(eq,
		Hill = c('-1 / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta*Beta / Alpha', '-1 / (1 + (Alpha / x)^Beta)^2*(Alpha / x)^Beta * log(Alpha / x)'),
		Hill_two = c('x / (Alpha + x)', '-Beta * x / (Alpha + x)^2'),
		Hill_three = c('-Gamma / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * Beta / Alpha', '-Gamma / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * log(Alpha / x)', '1 / (1 + (Alpha / x)^Beta)'),
		Hill_four = c('-(Gamma - Delta) / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * Beta / Alpha', '-(Gamma - Delta) / (1 + (Alpha / x)^Beta)^2 * (Alpha / x)^Beta * log(Alpha / x)', 
						'1 /(1 + (Alpha / x)^Beta)', '1 - 1 / (1 + (Alpha / x)^Beta)'),
		# Delta for Alpha_one and Epsilon for Beta_one
		Hill_five = c('(Gamma-1)/(1+(Alpha/x)^Beta)^2*(Alpha/x)^Beta*Beta/Alpha*(1-1/(1+(Delta/x)^Epsilon))', '(Gamma-1)/(1+(Alpha/x)^Beta)^2*(Alpha/x)^Beta*log(Alpha/x)*(1-1/(1+(Delta/x)^Epsilon))',
						'-1/(1+(Alpha/x)^Beta)*(1-1/(1+(Delta/x)^Epsilon))', '-(1+(Gamma-1)/(1+(Alpha/x)^Beta))/(1+(Delta/x)^Epsilon)^2*(Delta/x)^Epsilon*Epsilon/Delta',
						'-(1+(Gamma-1)/(1+(Alpha/x)^Beta))/(1+(Delta/x)^Epsilon)^2*(Delta/x)^Epsilon*log(Delta/x)'),
		
		
		Weibull = c('exp(Alpha + Beta * log(x) / log(10)) * exp(-exp(Alpha + Beta * log(x) / log(10)))', 
					'log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp(-exp(Alpha + Beta * log(x) / log(10)))'),
		
		Weibull_three = c('Gamma * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))',
							'Gamma * log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))', '1 - exp( -exp(Alpha + Beta * log(x) / log(10)))'),
		
		Weibull_four = c(' -(Delta - Gamma) * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))',
							' -(Delta - Gamma) * log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp( -exp(Alpha + Beta * log(x) / log(10)))', 
							'1 -exp( -exp(Alpha + Beta * log(x) / log(10)))', 'exp( -exp(Alpha + Beta * log(x) / log(10)))'),
		
		Logit = c('1 / (1 + exp(-Alpha - Beta * log(x) / log(10)))^2 * exp(-Alpha - Beta * log(x) / log(10))', 
					'1 / (1 + exp(-Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp(-Alpha - Beta * log(x) / log(10))'),
		
		Logit_three = c('Gamma / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * exp( -Alpha - Beta * log(x) / log(10))',
							'Gamma / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp( -Alpha - Beta * log(x) / log(10))',
							'1 / (1 + exp( -Alpha - Beta * log(x) / log(10)))'), 

		Logit_four = c('( -Delta + Gamma) / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * exp( -Alpha - Beta * log(x) / log(10))',
						'( -Delta + Gamma) / (1 + exp( -Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp( -Alpha - Beta * log(x) / log(10))',
						'1 / (1 + exp( -Alpha - Beta * log(x) / log(10)))', '1 - 1 / (1 + exp( -Alpha - Beta * log(x) / log(10)))'), 
		BCW = c('exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))', 
				'(x^Gamma - 1) / Gamma * exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))', 
				'(Beta * x^Gamma * log(x) / Gamma - Beta * (x^Gamma - 1) / Gamma^2) * exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))'),
		BCL = c('1 /(1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma))^2 * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)', 
				'1 / (1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma))^2 * (x^Gamma - 1) / Gamma * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)', 
				'-1 / (1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)) ^ 2 * (-Beta * x^Gamma * log(x) / Gamma + Beta * (x^Gamma - 1) / Gamma^2) * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)'),
		GL = c('1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * Gamma * exp(-Alpha - Beta * log(x) / log(10)) / (1 + exp(-Alpha - Beta * log(x) / log(10)))', 
				'1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * Gamma * log(x) / log(10) * exp(-Alpha - Beta * log(x) / log(10)) / (1 + exp(-Alpha - Beta * log(x) / log(10)))', 
				'-1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * log(1 + exp(-Alpha - Beta * log(x) / log(10)))'),
	
		Brain_Consens = c('-x / (1 + exp(Beta * Gamma) * x^Beta)', 
							'(1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)^2 * (Gamma * exp(Beta * Gamma) * x^Beta + exp(Beta * Gamma) * x^Beta * log(x))', 
							'(1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)^2 * Beta * exp(Beta * Gamma) * x^Beta'),
		
		BCV = c('-(1  +  Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)', 
				'-Alpha * x / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta) - 2 * Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)^2 * Gamma * (x / Gamma)^Delta',
				'Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)^2 * (2 * Beta * (x / Gamma)^Delta - (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta * Delta / Gamma)', 
				'Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)^2 * (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta * log(x / Gamma)'),
		  
		Cedegreen = c('-exp(-1 / (x^Beta)) / (1 + exp(Gamma * (log(x) - log(Delta))))',
						'-Alpha / (x^Beta) * log(x) * exp(-1 / (x^Beta)) / (1 + exp(Gamma * (log(x) - log(Delta))))',
						'(1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))^2 * (log(x) - log(Delta)) * exp(Gamma * (log(x) - log(Delta)))', 
						'-(1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))^2 * Gamma / Delta * exp(Gamma * (log(x) - log(Delta)))'),

		Beckon = c('(1 - 1 / (1 + (Beta / x)^Gamma)) / (1 + (x / Delta)^Epsilon)', 
					'Alpha / (1 + (Beta / x)^Gamma)^2 * (Beta / x)^Gamma * Gamma / Beta / (1 + (x / Delta)^Epsilon)',
					'Alpha / (1 + (Beta / x)^Gamma)^2 * (Beta / x)^Gamma * log(Beta / x) / (1 + (x / Delta)^Epsilon)', 
					'(Alpha + 1 - Alpha / (1 + (Beta / x)^Gamma)) / (1 + (x / Delta)^Epsilon)^2 * (x / Delta)^Epsilon * Epsilon / Delta',
					'-(Alpha + 1 - Alpha / (1 + (Beta / x)^Gamma)) / (1 + (x / Delta)^Epsilon)^2 * (x / Delta)^Epsilon * log(x / Delta)'),

		Biphasic = c('1 - 1 / (1 + 10^((x-Beta) * Gamma)) - 1 / (1 + 10^((Delta - x) * Epsilon))', 
						'-Alpha / (1 + 10^((x-Beta) * Gamma))^2 * 10^((x - Beta) * Gamma) * Gamma * log(10)',
						'Alpha / (1 + 10^((x - Beta) * Gamma))^2 * 10^((x - Beta) * Gamma) * (x - Beta) * log(10)',
						'-(1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))^2 * 10^((Delta - x) * Epsilon) * Epsilon * log(10)',
						'-(1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))^2 * 10^((Delta - x) * Epsilon) * (Delta - x) * log(10)')
	)
	
	for (i in seq(mpara)) jac[, i] <- eval(parse(text = jacFun[i]))
	return(jac)
}