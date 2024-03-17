#include "rayt.h"

namespace rayt {

	class Shape;
	class Material;
	typedef std::shared_ptr<Shape> ShapePtr;
	typedef std::shared_ptr<Material> MaterialPtr;

	//-------------------------------------------------------------------------

	class HitRec {
	public:
		float t;
		vec3 p;
		vec3 n;
		MaterialPtr mat;
	};

	//-------------------------------------------------------------------------

	class ScatterRec {
	public:
		Ray	ray;
		vec3 albedo;
	};

	class Material {
	public:
		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const = 0;
	};

	class Lambertian : public Material {
	public:
		Lambertian(const vec3& c)
			: m_albedo(c) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			vec3 target = hrec.p + hrec.n + random_in_unit_sphere();
			srec.ray = Ray(hrec.p, target - hrec.p);
			srec.albedo = m_albedo;
			return true;
		};

	private:
		vec3 m_albedo;
	};

	//-------------------------------------------------------------------------

	class Shape {
	public:
		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const = 0;
	};

	class ShapeList : public Shape {
	public:
		ShapeList() {}

		void add(const ShapePtr& shape) {
			m_list.push_back(shape);
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			HitRec temp_rec;
			bool hit_anything = false;
			float closest_so_far = t1;
			for (auto& p : m_list) {
				if (p->hit(r, t0, closest_so_far, temp_rec)) {
					hit_anything = true;
					closest_so_far = temp_rec.t;
					hrec = temp_rec;
				}
			}
			return hit_anything;
		}

	private:
		std::vector<ShapePtr> m_list;
	};

	class Sphere : public Shape {
	public:
		Sphere() {}
		Sphere(const vec3& c, float r, const MaterialPtr& mat)
			: m_center(c)
			, m_radius(r)
			, m_material(mat) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			vec3 oc = r.origin() - m_center;
			float a = dot(r.direction(), r.direction());
			float b = 2.0f*dot(oc, r.direction());
			float c = dot(oc, oc) - pow2(m_radius);
			float D = b*b - 4 * a*c;
			if (D > 0) {
				float root = sqrtf(D);
				float temp = (-b - root) / (2.0f*a);
				if (temp < t1 && temp > t0) {
					hrec.t = temp;
					hrec.p = r.at(hrec.t);
					hrec.n = (hrec.p - m_center) / m_radius;
					hrec.mat = m_material;
					return true;
				}
				temp = (-b + root) / (2.0f*a);
				if (temp < t1 && temp > t0) {
					hrec.t = temp;
					hrec.p = r.at(hrec.t);
					hrec.n = (hrec.p - m_center) / m_radius;
					hrec.mat = m_material;
					return true;
				}
			}

			return false;
		}

	private:
		vec3 m_center;
		float m_radius;
		MaterialPtr m_material;
	};

	class Triangle : public Shape {
	public:
		Triangle() {}
		Triangle(const vec3& v0, const vec3& v1, const vec3& v2, const MaterialPtr& m)
			: m_v0(v0)
			, m_v1(v1)
			, m_v2(v2)
			, m_material(m) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {

			vec3 e1 = m_v1 - m_v0;
			vec3 e2 = m_v2 - m_v0;

			vec3 alpha = cross(r.direction(), e2);
			float det = dot(e1, alpha);
			if (det > -EPSILON && det < EPSILON) {
				return false;
			}

			float invdet = 1.0f / det;
			vec3 oc = r.origin() - m_v0;

			// u が 0 <= u <= 1 を満たす。
			float u = dot(alpha, oc) * invdet;
			if (u < 0.0f || u > 1.0f) {
				return false;
			}

			vec3 beta = cross(oc, e1);

			// v が 0 <= v <= 1 かつ u + v <= 1 を満たす。
			float v = dot(r.direction(), beta) * invdet;
			if (v < 0.0f || u + v > 1.0f) {
				return false;
			}
			
			float t = dot(e2, beta) * invdet;
			if (t < t0 || t > t1) {
				return false;
			}

			hrec.t = t;
			hrec.mat = m_material;
			hrec.p = r.at(t);
			hrec.n = normalize(cross(e1, e2));
			return true;
		}

	private:
		vec3 m_v0, m_v1, m_v2;
		MaterialPtr m_material;
	};

	//-------------------------------------------------------------------------
	class Scene {
	public:
		Scene(int width, int height, int samples)
			: m_image(std::make_unique<Image>(width, height))
			, m_backColor(0.2f)
			, m_samples(samples) {
		}

		void build() {

			// カメラ

			// z正方向から見る
			vec3 lookfrom(0, 0, 8);
			vec3 lookat(0, 0, 0);
			vec3 vup(0, 1, 0);

			// y正方向から見る
			// vec3 lookfrom(0, 10, 0);
			// vec3 lookat(0, 0, 0);
			// vec3 vup(0, 0, 1);

			// 綺麗な角度
			// vec3 lookfrom(0, 10, 5);
			// vec3 lookat(0, 0, 0);
			// vec3 vup(0, 0, 1);
			float aspect = float(m_image->width()) / float(m_image->height());
			m_camera = std::make_unique<Camera>(lookfrom, lookat, vup, 30, aspect);

			ShapeList* world = new ShapeList();
			world->add(std::make_shared<Sphere>(vec3(0, 0, -100.05), 100, std::make_shared<Lambertian>(vec3(0.8f, 0.8f, 0.0f))));
			world->add(std::make_shared<Triangle>(vec3(-4, -1, 0.5), vec3(0, -1, 0.5), vec3(-3, 3, 0.5), std::make_shared<Lambertian>(vec3(0.8f, 0.0f, 0.0f))));
			world->add(std::make_shared<Triangle>(vec3(0, -2, 0.5), vec3(4, -2, 0.5), vec3(1, 2, 0.5), std::make_shared<Lambertian>(vec3(0.8f, 0.0f, 0.0f))));
			m_world.reset(world);
		}

		vec3 color(const rayt::Ray& r, const Shape* world) {
			HitRec hrec;
			if (world->hit(r, 0.001f, FLT_MAX, hrec)) {
				ScatterRec srec;
				if (hrec.mat->scatter(r, hrec, srec)) {
					return mulPerElem(srec.albedo, color(srec.ray, world));
				}
				else {
					return vec3(0);
				}
			}
			return backgroundSky(r.direction());
		}

		vec3 background(const vec3& d) const {
			return m_backColor;
		}

		vec3 backgroundSky(const vec3& d) const {
			vec3 v = normalize(d);
			float t = 0.5f * (v.getY() + 1.0f);
			return lerp(t, vec3(1), vec3(0.5f, 0.7f, 1.0f));
		}

		void render() {

			build();

			int nx = m_image->width();
			int ny = m_image->height();
#pragma omp parallel for schedule(dynamic, 1) num_threads(NUM_THREAD)
			for (int j = 0; j<ny; ++j) {
				std::cerr << "Rendering (y = " << j << ") " << (100.0 * j / (ny - 1)) << "%" << std::endl;
				for (int i = 0; i<nx; ++i) {
					vec3 c(0);
					for (int s = 0; s<m_samples; ++s) {
						float u = (float(i) + drand48()) / float(nx);
						float v = (float(j) + drand48()) / float(ny);
						Ray r = m_camera->getRay(u, v);
						c += color(r, m_world.get());
					}
					c /= m_samples;
					m_image->write(i, (ny - j - 1), c.getX(), c.getY(), c.getZ());
				}
			}

			stbi_write_bmp("render.bmp", nx, ny, sizeof(Image::rgb), m_image->pixels());
		}

	private:
		std::unique_ptr<Camera> m_camera;
		std::unique_ptr<Image> m_image;
		std::unique_ptr<Shape> m_world;
		vec3 m_backColor;
		int m_samples;
	};

} // namespace rayt

int main()
{
	int nx = 300;
	int ny = 150;
	int ns = 100;
	std::unique_ptr<rayt::Scene> scene(std::make_unique<rayt::Scene>(nx, ny, ns));
	scene->render();
	return 0;
}

