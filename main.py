import wx
from blockmodel import *
from numpy import *
from density import *

#used to display the result of running a block modeling algorithm,
#including the agents in each block, T-value, and matrices displaying the 
#density of ties between blocks for each tie matrix in the stack
class Display(wx.Frame):
    def __init__(self, title):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=title, size=(650,600)) 
        
        self.text = wx.TextCtrl(self)

#the main app window
class Frame(wx.Frame):
    def __init__(self, title):
        wx.Frame.__init__(self, None, title=title, size=(700,150))        
        self.Center()

        self.main_panel = wx.Panel(self)

        #panel with buttons and text boxes
        self.info_panel = wx.Panel(self.main_panel)
        self.info_sizer = wx.FlexGridSizer(rows=0,cols=2)
        self.info_sizer.AddGrowableCol(1, proportion=5)
        self.info_panel.SetSizer(self.info_sizer)

        #buttons panel
        button_panel = wx.Panel(self.main_panel, pos = (0, 20))
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        button_panel.SetSizer(button_sizer)

        #buttons
        self.concor = wx.Button(button_panel, -1, "CONCOR")
        self.random1 = wx.Button(button_panel, -1, "Best of 100 (before optimizing)")
        self.random2 = wx.Button(button_panel, -1, "Best of 100 (optimizing each time")
        self.optimize = wx.Button(button_panel, -1, "Optimize")
        self.clear = wx.Button(button_panel, -1, "Clear Blocks")

        self.concor.Bind(wx.EVT_BUTTON, self.on_concor)
        self.random1.Bind(wx.EVT_BUTTON, self.on_random1)
        self.random2.Bind(wx.EVT_BUTTON, self.on_random2)
        self.optimize.Bind(wx.EVT_BUTTON, self.on_optimize)
        self.clear.Bind(wx.EVT_BUTTON, self.on_clear)

        button_sizer.AddMany([(self.concor, 0),
        						(self.random1, 0),
        						(self.random2, 0),
        						(self.optimize, 0),
        						(self.clear, 0)])

        #text boxes and parameter entry
        self.num_blocks = self.add_infobox("How many blocks do you want?")
        self.num_mats = self.add_infobox("How many matrices are in your stack?")
        self.enter = wx.Button(self.info_panel, -1, "Enter Matrix Data")
        self.info_sizer.Add(self.enter, flag = wx.ALL, border = 5)
        self.enter.Bind(wx.EVT_BUTTON, self.matrix_entry)

        #the blockmodel display frame
        self.view = Display('Block Model')

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(button_panel, flag = wx.BOTTOM, border = 5)
        main_sizer.Add(self.info_panel, 0, wx.EXPAND)
        main_sizer.Add(self.view, 1, wx.EXPAND)
        self.main_panel.SetSizer(main_sizer)

        self.model = block_mat([])

    #for creating parameter boxes
    def add_infobox(self, title):
        infolabel = wx.StaticText(self.info_panel, label=title)
        infobox = wx.TextCtrl(self.info_panel)
        self.info_sizer.Add(infolabel, -1, wx.ALIGN_RIGHT|wx.ALL, 5)
        self.info_sizer.Add(infobox, -1, wx.EXPAND|wx.ALL, 5)
        return infobox

    #for getting matrix data from the user
    def matrix_entry(self, event):
        n = self.num_mats.GetValue()
        n2 = self.num_blocks.GetValue()

        try:
            n2 = int(n2)
            self.model.num_blocks = n2
        except:
            self.num_blocks.SetValue("Please enter a positive integer input.")

        try:
            n = float(n)
            m = int(n)
            self.entry = DataDialog('Enter a Matrix', m)
            if n != m or n <= 0:
                self.num_mats.SetValue("Enter a positive integer, please!")
            else:
                self.entry.Show()
        except:
            self.num_mats.SetValue("Enter a positive integer, please!")

    #find a blocking based on the CONCOR algorithm
    def on_concor(self, event):
        self.model.concor()

        text = "Blocks from CONCOR: \n"
        for i in range(self.model.num_blocks):
            n = len(self.model.blocks[i])
            text += "\n\n Block " + str(i) + ": "
            for p in self.model.blocks[i][0: n - 1]:
                text += str(p + 1) + ", "
            text += str(self.model.blocks[i][n - 1])
        text +=  "\n\n T-value: "  + str(self.model.T)
        text += self.display_density()

        self.view.text.SetValue(text)
        if not self.view.IsShown():
            self.view.Show()

    #find a blocking based from the best 100 randomly generated block models
    def on_random1(self, event):
        self.model.random_blocks()

        text =  "Best of 100 Random Blockings (Non-Optimized): \n" 
        for i in range(self.model.num_blocks):
            n = len(self.model.blocks[i])
            text +=  "\n\n Block " + str(i) + ": " 
            for p in self.model.blocks[i][0: n - 1]:
                text += str(p + 1) + ", "
            text += str(self.model.blocks[i][n - 1])
        text +=  "\n\n T-value: "  + str(self.model.T)
        text += self.display_density()

        self.view.text.SetValue(text)
        if not self.view.IsShown():
            self.view.Show()

    #run 100 trials of 1. generating a random block model and 2. optimizing it
    #return the optimal model
    def on_random2(self, event):
        first_random = self.model.single_random_block()
        self.model.blocking(first_random)
        self.model.optimize_by_parts(partial_T, T)

        old_T = self.model.T

        new_model = block_mat(self.model.matrix_list, self.model.num_blocks)

        for i in range(100):
            new_blocks = new_model.single_random_block()
            new_model.blocking(new_blocks)
            new_model.optimize_by_parts(partial_T, T)
            new_T = new_model.T

            if new_T > old_T:
                self.model.blocking(new_blocks)
                old_T = new_T
                self.model.T = new_T

        text =  "Best of 100 Randomly Generated and Optimized Blockings: \n" 
        for i in range(self.model.num_blocks):
            n = len(self.model.blocks[i])
            text +=  "\n\n Block " + str(i) + ": " 
            for p in self.model.blocks[i][0: n - 1]:
                text += str(p + 1) + ", "
            text += str(self.model.blocks[i][n - 1])
        text +=  "\n\n T-value: "  + str(self.model.T)
        text += self.display_density()

        self.view.text.SetValue(text)
        if not self.view.IsShown():
            self.view.Show()

    #optimize the current block model associated with the data
    def on_optimize(self, event):
        self.model.optimize()

        text =  "Optimized Blocks: \n" 
        for i in range(self.model.num_blocks):
            n = len(self.model.blocks[i])
            text +=  "\n\n Block " + str(i) + ": " 
            for p in self.model.blocks[i][0: n - 1]:
                text += str(p + 1) + ", "
            text += str(self.model.blocks[i][n - 1])
        text +=  "\n\n T-value: "  + str(self.model.T)
        text += self.display_density()

        self.view.text.SetValue(text)
        if not self.view.IsShown():
            self.view.Show()

    #clear data and block model
    def on_clear(self, event):
        self.view.Destroy()
        self.model.num_blocks = 1
        self.num_blocks.SetValue("")
        self.num_mats.SetValue("")
        self.model.blocking([[]])

    #used for displaying the matrices corresponding to the density of ties between blocks
    def display_density(self):
        text = ""
        block_mats = blocked_mats(self.model.blocks, self.model.matrix_list, self.model.ntwk_size)
        for i in range(self.model.num_mats):
            text +=  "\n\n Density Matrix (" + self.model.matrix_names[i] + "): " 
            density_mat = block_dense(self.model.blocks, block_mats[i])
            for row in density_mat:
                text += "\n" + str(row)

        return text

#used for gett matrix data from the user
class DataDialog(wx.Frame):
    def __init__(self, title, number_left):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=title, size=(700,200))

    	self.main_panel = wx.Panel(self)

        self.info_panel = wx.Panel(self.main_panel)
        self.info_sizer = wx.FlexGridSizer(rows=0,cols=2)
        self.info_sizer.AddGrowableCol(1, proportion=5)
        self.info_panel.SetSizer(self.info_sizer)

    	self.left = number_left
    	self.done = 0

        self.info = wx.StaticText(self.info_panel, label = "Enter title of Matrix # " + str(self.done + 1) + ":")
    	self.labelbox = wx.TextCtrl(self.info_panel)
    	self.info_sizer.Add(self.info, -1, wx.ALIGN_LEFT|wx.ALL, 5)
    	self.info_sizer.Add(self.labelbox, -1, wx.EXPAND|wx.ALL, 5)

    	#add multiline textbox for matrix entry
    	self.matrix_info = wx.StaticText(self.info_panel, -1, "Paste or type matrix with spaces separating elements in each row, and place each row on a new line.")
        self.matrix_info.Wrap(200)
    	self.matrix_text = wx.TextCtrl(self.info_panel, -1, "Example: \n1 0\n0 1 \nwould denote a 2x2 identity matrix", size = (100,100), style = wx.TE_MULTILINE)
    	self.matrix_text.SetInsertionPoint(0)
    	self.info_sizer.Add(self.matrix_info, -1, wx.ALIGN_LEFT|wx.ALL, 5)
    	self.info_sizer.Add(self.matrix_text, -1, wx.EXPAND|wx.ALL, 5)

    	self.enterbutton = wx.Button(self.info_panel, -1, "Set Data")
    	self.info_sizer.Add(self.enterbutton, -1, wx.ALL, 5)
        self.enterbutton.Bind(wx.EVT_BUTTON, self.get_input)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.info_panel, 0, wx.EXPAND)
        self.main_panel.SetSizer(main_sizer)

    def get_matrix(self):
    	string_version = str(self.matrix_text.GetValue())
        list_version = string_version.split('\n')

        for i in range(len(list_version)):
            list_version[i] = list_version[i].split(' ')
            for j in range(len(list_version[i])):
                el = list_version[i][j]

                try: el = float(el)
                except: self.matrix_text.SetValue("Please enter matrix data as instructed.")

                if int(el) != el:
                    self.matrix_text.SetValue("Please enter matrix data as instructed.")
                else: list_version[i][j] = int(el)

        n = len(list_version)

    	if n != len(list_version[0]):
    		self.matrix_text.SetValue("This matrix isn't square!")
        else:
            for row in list_version:
                if len(row) != n:
                    self.matrix_text.SetValue("The rows aren't all the same size!")
                else: return list_version

    def get_title(self):
    	return self.labelbox.GetValue()

    def get_input(self, event):
        parent = self.GetParent()
        
        if not self.get_title():
            self.labelbox.SetValue("Please enter a title for the matrix as instructed.")
        else:
            mat = self.get_matrix()
            ones_mat = find_binary(mat, .15, .20)
            parent.model.matrix_list.append(ones_mat)
            parent.model.matrix_names.append(self.get_title())

            if self.left == 1:
                parent.model.update_info()
                self.Destroy()
            else:
                self.done += 1
                self.left -= 1
                self.info.SetLabel("Enter title of Matrix # " + str(self.done + 1) + ":")
                self.info_sizer.Layout()
                self.matrix_text.SetValue('Next matrix,please!')
                self.labelbox.SetValue("")

app = wx.App()

frame = Frame('Block Model')
frame.Show()

app.MainLoop()
