import streamlit as st

st.set_page_config(layout="wide", page_title='DNA sequencying')

import streamlit.components.v1 as htmlviewer
# Title Msg#1
st.title('DNA로부터 단백질까지')

with open('./index.html','r', encoding='utf-8') as f:
    html = f.read()
with open('./index2.html', 'r', encoding='utf-8') as f2:
    html2 = f2.read()
    f.close()


# html = '''
# <html>
#     <head>
#         <title> This is Jihye's html </title>
#     </head>

#     <body>
#         <h1>Topic</h1>
#         <h2>SubTopic</h2>
#     </body>
# </html>
# '''

# Box#1(4), Box#2(1)
col1, col2 = st.columns((6,1))
with col1:
    with st.expander('Contents #1'):
        url = 'https://youtu.be/89YbXYPu50c?si=ucF84JOB_Gsn_vRd'
        st.info('참고 영상: DNA란 무엇인가?')
        st.video(url)

    with st.expander('Contents #2'):
        #st.write(html, unsafe_allow_html=True)
        htmlviewer.html(html, height=800, width=1500, scrolling=True)

    with st.expander('Contents #3'):
        #st.write(html, unsafe_allow_html=True)
        htmlviewer.html(html2, height=1200, width=1500, scrolling=True)


with col2:
    with st.expander('Tips..'):
        st.info('Tips..')
st.markdown('<hr>', unsafe_allow_html=True)
st.write('<font color="BLUE">(c)copyright. all rights reserved by Jihye Kim', unsafe_allow_html=True)