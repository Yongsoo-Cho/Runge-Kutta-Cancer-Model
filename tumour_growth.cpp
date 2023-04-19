float tumour_growth(float tg, float rad, float delta, float tau, float gamma)
{
	return (delta*rad) - (tg/tau) - (gamma* (tg*tg));
}
