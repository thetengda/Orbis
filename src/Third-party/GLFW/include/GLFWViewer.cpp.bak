#include "GLFWViewer.h"
#include <iostream>

void MyPerspective(GLdouble fov, GLdouble aspectRatio, GLdouble zNear, GLdouble zFar);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
 //map<GLFWwindow*, GLFWViewer*> GLFWViewer::mapptr;
 //GLFWViewer* getViewerPtr(GLFWwindow* window);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
// Different Objects
void MyCylinder(GLdouble r, GLdouble l, int edgenum);
void MyFrame();
void MyAerial(GLfloat size);

vector<GLFWViewer*> GLFWViewer::vptr;
GLFWViewer viewer ;

 using namespace std;

 bool mb_OK;
 bool mb_MPs = true;
 bool mb_trajetory = true;
 bool mb_KFs = true;
 bool mb_overlook = false;
 bool mb_grid = true;
 bool mb_axis = true;
 double m_interval = 0.1;
 double m_dist = 5;
 double m_scale = 1;
 double m_rot_X = 0;
 double m_rot_Y = 0;
 double m_rot_Z = 0;
 double m_rot_X0 = 0;
 double m_rot_Y0 = 0;
 double m_rot_Z0 = 0;
 double m_trans_X = 0;
 double m_trans_Y = 0;
 double m_trans_X0 = 0;
 double m_trans_Y0 = 0;

 double m_mouseL_x = 0;
 double m_mouseL_y = 0;
 double m_mouseR_x = 0;
 double m_mouseR_y = 0;
 double m_dx = 0;
 double m_dy = 0;
 bool is_pressL = false;
 bool is_pressR = false;
 bool mb_stop=false;

GLFWViewer::GLFWViewer()
{
	GLFWViewer::vptr.push_back(this);
}


GLFWViewer::~GLFWViewer()
{
	GLFWViewer::vptr.erase(std::find(GLFWViewer::vptr.begin(), GLFWViewer::vptr.end()-1,this));

}



//GLFWViewer* getViewerPtr(GLFWwindow* window)
//{
//	for (auto ite = GLFWViewer::vptr.begin(); ite != GLFWViewer::vptr.end(); ite++)
//	{
//		if ((*ite)->window == window)
//			return *ite;
//	}
//	return nullptr;
//}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (!mb_overlook)
		MyPerspective(45, (float)width / height, 0.01, 100.0);
	else glOrtho(-3.0*width / height, 3.0*width / height, -3, 3, -10, 1000);
}
void GLFWViewer::Run()
{
	/* Initialize the library */
	if (!glfwInit())
		return;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(100, 100, "Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return;
	}
	/* Make the window's context current */
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetKeyCallback(window, key_callback);
	glfwSetWindowSize(window,800 , 600);
	glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, 1);
	
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	 //   glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHTING);
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window)&&!mb_stop)
	{
		Sleep(20);
		unique_lock<mutex> lock(m_mutex);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//paint
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(m_trans_X, m_trans_Y, -m_dist * sqrt(3));

		if (!mb_overlook)
		{
			//gluLookAt(m_dist,m_dist,m_dist,0,0,0,0,0,1);
			glRotatef(-50 + m_rot_X, 1, 0, 0);
			glRotatef(-135, 0, 0, 1);
			glRotatef(0, 1, 0, 0);
			glRotatef(m_rot_Y, 0, 1, 0);
			glRotatef(m_rot_Z, 0, 0, 1);
		}
		else
		{
			//gluLookAt(0,0,m_dist,0,0,0,0,1,0);
		}

		m_interval = 0.1;

		glScalef(m_scale, m_scale, m_scale);
		double temp = 10 / m_scale;
		while (temp > 5)
		{
			temp /= 5;
			m_interval *= 5;
			if (temp > 2)
			{
				temp /= 2;
				m_interval *= 2;
			}
		}
		if (mb_grid)
		{
			for (double  i = 0; i < 100 * (1 / m_scale); i += m_interval)
			{
				double z = 1.5*sqrt(sqrt(0.00002*i / m_interval));
				if (i == 0) z = 1.5*sqrt(sqrt(0.00002 * 1));
				glColor3f(0.70f + z, 0.70f + z, 0.70f + z);
				glLineWidth(1.0);
				glBegin(GL_LINES);
				glVertex3f(-100 * (1 / m_scale), i, 0);
				glVertex3f(100 * (1 / m_scale), i, 0);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f(-100 * (1 / m_scale), -i, 0);
				glVertex3f(100 * (1 / m_scale), -i, 0);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f(i, -100 * (1 / m_scale), 0);
				glVertex3f(i, 100 * (1 / m_scale), 0);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f(-i, -100 * (1 / m_scale), 0);
				glVertex3f(-i, 100 * (1 / m_scale), 0);
				glEnd();
			}
		}
		if (mb_axis)
		{
			glPushMatrix();
			glColor3f(0.0f, 0.0f, 1.0f);
			MyCylinder(0.025f * 1 / m_scale, 1.5*m_interval, 8);
			glPopMatrix();
			glPushMatrix();
			glRotatef(-90, 1, 0, 0);
			glColor3f(0.0f, 1.0f, 0.0f);
			MyCylinder(0.025f * 1 / m_scale, 1.5*m_interval, 8);
			glPopMatrix();
			glPushMatrix();
			glRotatef(90, 0, 1, 0);
			glColor3f(1.0f, 0.0f, 0.0f);
			MyCylinder(0.025f * 1 / m_scale, 1.5*m_interval, 8);
			glPopMatrix();
		}
		

		if (mb_trajetory)
		{
			float c_s[3] = { 250/255.0,167/ 255.0,85/ 255.0 };
			float c_e[3] = { 80/255.0,183/255.0,193/255.0 };
			glColor4f(c_s[0],c_s[1],c_s[2],0.5);
			glLineWidth(6.0f);
			for (const auto &t : mv_trajectory)
			{
				glBegin(GL_LINE_STRIP);
				for ( int i=0;i<t.size();i++)
				{
					auto p = t[i];
					int n = 20-(t.size()-1 - i);
					if (i > t.size() - 20 && i <=t.size()-5)
					{
						glColor4f(
							c_s[0] + (c_e[0] - c_s[0]) / 15.0 * n,
							c_s[1] + (c_e[1] - c_s[1]) / 15.0 * n,
							c_s[2] + (c_e[2] - c_s[2]) / 15.0 * n,0.5);
					}
					glVertex3d(p(0), p(1), p(2));
				}
				glEnd();
			}

		}
		if (mb_KFs)
		{
			for (const auto &fs : mv_frames)
			{
				for (const auto &f : fs)
				{
					Eigen::Matrix4d Twb = Eigen::Matrix4d::Identity();
					Twb.block<3, 3>(0, 0) = f.first;
					Twb.block<3, 1>(0, 3) = f.second;
					glPushMatrix();

					glMultMatrixd(Twb.data());
					
					MyAerial(0.1);

					glPopMatrix();
				}
			}
		}
		if (mb_MPs)
		{
			double c_l[3] = {0,0,0.2};
			double c_h[3] = { 0.9,0.9,1 };
			double low = -5, high = 5;
			glPointSize(4.0f);
			for (const auto &pc : mv_pointCloud)
			{
				for (const auto &p:pc)
				{
					double h = (p(2) - (-5))/10;
					glColor3f(c_l[0]+(c_h[0]-c_l[0])*h,
						c_l[1] + (c_h[1] - c_l[1])*h, 
						c_l[2] + (c_h[2] - c_l[2])*h);
					glBegin(GL_POINTS);
					glVertex3f(p(0),p(1),p(2));
					glEnd();
				}
			}
			
		}
		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();
	//cout<<window->
	//delete window;
}

void GLFWViewer::ClearView()
{
	unique_lock<mutex> lock(m_mutex);
	mv_trajectory.clear();
	mv_pointCloud.clear();
}

void GLFWViewer::SetFrames(const Frames &f)
{
	unique_lock<mutex> lock(m_mutex);
	mv_frames.clear();
	mv_frames.push_back(f);
}

void GLFWViewer::SetPointCloud(const PointCloud &pc)
{
	unique_lock<mutex> lock(m_mutex);
	mv_pointCloud.clear();
	mv_pointCloud.push_back(pc);
}

void GLFWViewer::SetFrames(const vector<Frames> &vf)
{
	unique_lock<mutex> lock(m_mutex);
	mv_frames = vf;
}

void GLFWViewer::SetPointCloud(const vector<PointCloud> &vpc)
{
	unique_lock<mutex> lock(m_mutex);
	mv_pointCloud = vpc;
}

void GLFWViewer::SetTrajectory(const Trajectory & t)
{
	unique_lock<mutex> lock(m_mutex);
	mv_trajectory.clear();
	mv_trajectory.push_back(t);
}

void GLFWViewer::SetTrajectory(const vector<Trajectory>& vt)
{
	unique_lock<mutex> lock(m_mutex);
	mv_trajectory = vt;
}

void GLFWViewer::Show()
{
	if(t) return;
	mb_stop = false;
	t = new thread(&GLFWViewer::Run, this);
	//t->detach();
}

void GLFWViewer::Hide()
{
	if(!t) return;

	mb_stop = true;
	//Sleep(2000);
	t->join();
	delete t;
	window = nullptr;
	t = nullptr;
}

void MyCylinder(GLdouble r, GLdouble l, int edgenum)
{
	for (int i = 0; i < edgenum; i++)
	{
		double a1 = i * 2 * 3.1416 / edgenum;
		double a2 = (i + 1) * 2 * 3.1416 / edgenum;

		double x1 = r * cos(a1);
		double y1 = r * sin(a1);
		double x2 = r * cos(a2);
		double y2 = r * sin(a2);
		glBegin(GL_TRIANGLES);
		glVertex3f(x1, y1, 0);
		glVertex3f(x2, y2, 0);
		glVertex3f(x1, y1, l);
		glEnd();

		glBegin(GL_TRIANGLES);
		glVertex3f(x1, y1, l);
		glVertex3f(x2, y2, l);
		glVertex3f(x2, y2, 0);
		glEnd();
	}
}

void MyPerspective(GLdouble fov, GLdouble aspectRatio, GLdouble zNear, GLdouble zFar)
{
	// 使用glu库函数，需要添加glu.h头文件
	//gluPerspective( fov, aspectRatio, zNear, zFar );
	// 使用OpenGL函数，但是需要添加math.h头文件
	GLdouble rFov = fov * 3.14159265 / 180.0;
	glFrustum(-zNear * tan(rFov / 2.0) * aspectRatio,
		zNear * tan(rFov / 2.0) * aspectRatio,
		-zNear * tan(rFov / 2.0),
		zNear * tan(rFov / 2.0),
		zNear, zFar);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		m_mouseL_x = xpos;
		m_mouseL_y = ypos;
		m_rot_Z0 = m_rot_Z;
		m_rot_X0 = m_rot_X;
		is_pressL = true;
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
	{
		m_mouseR_x = 0;
		m_mouseR_y = 0;
		is_pressL = false;
	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		m_mouseR_x = xpos;
		m_mouseR_y = ypos;
		m_trans_X0 = m_trans_X;
		m_trans_Y0 = m_trans_Y;
		is_pressR = true;
	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
	{
		m_mouseR_x = 0;
		m_mouseR_y = 0;
		is_pressR = false;
	}
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (is_pressL)
	{
		m_rot_Z = m_rot_Z0 + (xpos - m_mouseL_x)*0.2;
		m_rot_X = m_rot_X0 + (ypos - m_mouseL_y)*0.2;
	}
	if (is_pressR)
	{
		m_trans_Y = m_trans_Y0 - (ypos - m_mouseR_y)*0.015;
		m_trans_X = m_trans_X0 + (xpos - m_mouseR_x)*0.015;
	}
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	if (yoffset > 0) m_scale *= 1.1;
	if (yoffset < 0) m_scale /= 1.1;
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_O && action == GLFW_PRESS)
		mb_overlook = (bool)(1-(int)mb_overlook);
}

void MyAerial(GLfloat size)
{
	Eigen::Matrix4d T;
	T << 0, 0, 1, 0,
		0, 1, 0, 0,
		1, 0, 0, 0,
		0, 0, 0, 1;
	glPushMatrix();

	glMultMatrixd(T.data());

	float h = size/2;
	float l = sqrt( 2)*size;
	float r = size*0.7 ;
	glColor3f(0.2,0.2,0.8);
	glLineWidth(4.0f);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(size, size, h);
	glVertex3d(0, 0, 0);
	glVertex3d(-size, -size, h);
	glVertex3d(0, 0, 0);
	glVertex3d(-size, size, h);
	glVertex3d(0, 0, 0);
	glVertex3d(size, -size, h);
	glEnd();
	for (int i = 0; i < 4; i++)
	{
		double dd1 = 90.0/180*3.1415926;
		int n = 20;
		double dd2 = 360.0 / n / 180 * 3.1415926;
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j <= n; j++)
		{
			glVertex3d(cos(dd1*(i+0.5))*l+cos(dd2*j)*r,sin(dd1*(i + 0.5))*l+sin(dd2*j)*r,h);
		}
		glEnd();
	}
	glPopMatrix();

}

void MyFrame()
{
	float w = 0.1;
	float h = 0.1;
	float z = w * 0.6;
	glLineWidth(2);
	glColor3f(0.3f, 0.3f, 1.0f);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(w, h, z);
	glVertex3f(0, 0, 0);
	glVertex3f(w, -h, z);
	glVertex3f(0, 0, 0);
	glVertex3f(-w, -h, z);
	glVertex3f(0, 0, 0);
	glVertex3f(-w, h, z);

	glVertex3f(w, h, z);
	glVertex3f(w, -h, z);

	glVertex3f(-w, h, z);
	glVertex3f(-w, -h, z);

	glVertex3f(-w, h, z);
	glVertex3f(w, h, z);

	glVertex3f(-w, -h, z);
	glVertex3f(w, -h, z);
	glEnd();
}